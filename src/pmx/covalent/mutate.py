"""Emit a pmx hybrid residue using the graph-MCS mapper.

This wires :func:`pmx.covalent.map_residues` (graph-guaranteed mapping) into the
pmx dual-topology hybrid emit, reusing the proven output machinery from
:mod:`pmx.nsaa` (dummy generation, bonded-term merging, RTP/MTP writers).  It is
the emit-path counterpart of :func:`pmx.nsaa.mutate_nsaa`, differing only in the
mapping step -- so the resulting hybrid maps atoms by molecular topology rather
than by a 3-D MCS plus the polar-H patch, and common atoms such as ``HN`` are
mapped by identity instead of being dropped to dummies.

Two entry points:

* :func:`build_hybrid_residue` -- residue-level: given two *parameterized*
  residues (atom types/charges/masses + bond connectivity), returns the hybrid
  RTP/MTP text.  ``TargetResidue.from_rtp`` produces a suitable target.
* :func:`mutate` -- model-level: fetch the WT residue from a protein Model,
  parameterize it from the force field, build the target, and write the hybrid
  structure + RTP/MTP.
"""

from __future__ import annotations

import os
from copy import deepcopy
from typing import Optional, Tuple

import numpy as np

from .. import library
from ..model import Model
from ..ffparser import RTPParser
from ..utils import get_ff_path, list2file
from ..library import pmx_data_file

from .mapping import map_residues, attach_bonds_from_rtp, _element_of
from .target import TargetResidue

# Reuse the battle-tested emit helpers from nsaa (import, do not modify).
from ..nsaa import (
    _merge_molecules, _make_bstate_dummies, _make_transition_dics,
    _update_bond_lists, _improps_as_atoms, _generate_dihedral_entries,
    _generate_improp_entries, _make_rtp, _make_mtp, _write_atp_fnb,
    _rtp_terminal, _assign_rtp_entries, _assign_mass_atp,
    _validate_atom_parameters, _get_residue_types, _ARN_CONTENT,
)


# ── coordinate placement ──────────────────────────────────────────────────────

def _kabsch(P: np.ndarray, Q: np.ndarray):
    """Rigid transform (R, t) minimizing ||R P + t - Q|| (Kabsch)."""
    Pc = P - P.mean(0)
    Qc = Q - Q.mean(0)
    H = Pc.T @ Qc
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    t = Q.mean(0) - R @ P.mean(0)
    return R, t


def _superpose_target(residue1, residue2, pairs) -> None:
    """Superpose *residue2* onto *residue1* using mapped heavy-atom pairs.

    Mapped atoms end up (near-)coincident with their WT partners; unmapped
    target atoms (dummies) are carried into the WT frame so they get sensible
    coordinates in the emitted structure.
    """
    heavy = [(i, j) for i, j in pairs
             if _element_of(residue1.atoms[i]) != "H"
             and _element_of(residue2.atoms[j]) != "H"]
    if len(heavy) < 3:
        heavy = pairs  # fall back to all pairs if too few heavy anchors
    if len(heavy) < 3:
        return
    Q = np.array([residue1.atoms[i].x for i, _ in heavy], dtype=float)
    P = np.array([residue2.atoms[j].x for _, j in heavy], dtype=float)
    R, t = _kabsch(P, Q)
    for a in residue2.atoms:
        a.x = list(R @ np.array(a.x, dtype=float) + t)


# ── force-field helpers ───────────────────────────────────────────────────────

def _load_rtp(ff: str):
    ffpath = get_ff_path(ff)
    rtp_file = os.path.join(ffpath, "merged.rtp")
    if not os.path.isfile(rtp_file):
        rtp_file = os.path.join(ffpath, "aminoacids.rtp")
    return RTPParser(rtp_file), ffpath


def _inter_residue_bonds(rtp, resname, term):
    entr = _rtp_terminal(rtp[resname], term) if term else rtp[resname]
    return [[a1, a2] for a1, a2 in entr["bonds"]
            if a1.startswith(("+", "-")) or a2.startswith(("+", "-"))]


# ── residue-level hybrid build ────────────────────────────────────────────────

def build_hybrid_residue(residue1, residue2, ff: str, mut_resname: str,
                         term: Optional[str] = None,
                         write_ff_types: bool = True) -> Tuple[str, str, str]:
    """Build a hybrid residue from two parameterized residues via graph mapping.

    Both residues must already carry atom types/charges/masses and bond
    connectivity (``TargetResidue.from_rtp`` provides this for the target; a WT
    residue can be prepared with :func:`_assign_rtp_entries` +
    :func:`attach_bonds_from_rtp`).

    Returns ``(rtp_text, mtp_text, hybrid_name)``.  If *write_ff_types*, new
    ``DUM_`` atom types are appended to the force field's ``atomtypes.atp`` /
    ``ffnonbonded.itp``.
    """
    rtp, ffpath = _load_rtp(ff)
    bCharmm = "charmm" in ff

    _validate_atom_parameters(residue1, "state A")
    _validate_atom_parameters(residue2, "state B")
    _assign_mass_atp(residue1, residue2,
                     os.path.join(ffpath, "atomtypes.atp"), itp_file=None)

    # --- graph-MCS mapping (replaces domapping + polar-H patch) ---
    mapping = map_residues(residue1, residue2)
    pairs1, pairs2 = mapping.pairs1(), mapping.pairs2()

    # place the target (and its dummies) in the WT frame
    _superpose_target(residue1, residue2, list(zip(pairs1, pairs2)))

    # --- build hybrid residue (mirrors nsaa.mutate_nsaa) ---
    merged_atoms2 = []
    residue1.batoms = []
    atom_pairs = []
    for id1, id2 in zip(pairs1, pairs2):
        at1, at2 = residue1.atoms[id1], residue2.atoms[id2]
        at1.atomtypeB, at1.qB, at1.mB, at1.nameB = \
            at2.atomtype, at2.q, at2.m, at2.name
        residue1.batoms.append(at2)
        merged_atoms2.append(at2)
        atom_pairs.append([at1, at2])

    dummies = [a for a in residue2.atoms if a not in merged_atoms2]
    _merge_molecules(residue1, dummies)
    _make_bstate_dummies(residue1)

    if write_ff_types:
        _write_atp_fnb("atomtypes.atp", "ffnonbonded.itp", residue1, ffpath)

    _, badic = _make_transition_dics(atom_pairs, residue1)
    _update_bond_lists(residue1, badic)

    entr1 = _rtp_terminal(rtp[residue1.resname], term) if term \
        else rtp[residue1.resname]
    if residue2.resname in rtp:
        entr2 = _rtp_terminal(rtp[residue2.resname], term) if term \
            else rtp[residue2.resname]
    else:
        entr2 = {"atoms": [], "bonds": [], "diheds": [], "improps": [], "cmap": []}

    dih1 = _improps_as_atoms(entr1["diheds"], residue1, residue2)
    dih2 = _improps_as_atoms(entr2["diheds"], residue1, residue2, use_b=True)
    im1 = _improps_as_atoms(entr1["improps"], residue1, residue2)
    im2 = _improps_as_atoms(entr2["improps"], residue1, residue2, use_b=True)

    cmap = rtp[residue1.resname].get("cmap", []) if bCharmm else []
    dihi_list = _generate_dihedral_entries(dih1, dih2, residue1, atom_pairs)
    ii_list = _generate_improp_entries(im1, im2, residue1)

    bond_neigh = _inter_residue_bonds(rtp, residue1.resnA if hasattr(
        residue1, "resnA") else residue1.resname, term)

    one_letter = library._one_letter.get(residue1.resname, residue1.resname[0])
    hybrid_name = (one_letter + "2" + mut_resname.upper())[:4]
    residue1.resnA, residue1.resnB = residue1.resname, mut_resname
    residue1.set_resname(hybrid_name)

    rtp_text = "".join(_make_rtp(residue1, ii_list, dihi_list,
                                 bond_neigh, cmap, ff))
    mtp_text = "".join(_make_mtp(residue1, ii_list, dihi_list, rotations=None))
    return rtp_text, mtp_text, hybrid_name


# ── model-level mutate ────────────────────────────────────────────────────────

def mutate(model, mut_resid: int, target: str, ff: str,
           mut_chain: Optional[str] = None, inplace: bool = False,
           write_files: bool = True, verbose: bool = False):
    """Mutate ``mut_resid`` in *model* to an RTP-named *target* residue.

    Builds the hybrid structure and writes ``<hybrid>.rtp`` (into the force
    field dir) and ``<hybrid>.mtp`` (into the working dir), mirroring
    :func:`pmx.nsaa.mutate_nsaa` but using the graph-MCS mapper.

    Returns the hybrid residue name.
    """
    rtp, ffpath = _load_rtp(ff)
    bCharmm = "charmm" in ff

    m = model if inplace else deepcopy(model)
    if m.unity == "nm":
        m.nm2a()

    residue1 = m.fetch_residue(idx=mut_resid, chain=mut_chain)
    chain = residue1.chain
    term = ("C" if residue1 == chain.cterminus()
            else ("N" if residue1 == chain.nterminus() else None))

    # parameterize + bond the WT residue from the force field
    _assign_rtp_entries(residue1, rtp, term)
    if not any(getattr(a, "bonds", None) for a in residue1.atoms):
        attach_bonds_from_rtp(residue1, ff)

    # build the target residue from its RTP entry
    residue2 = TargetResidue.from_rtp(target, ff).residue

    rtp_text, mtp_text, hybrid_name = build_hybrid_residue(
        residue1, residue2, ff, target, term=term, write_ff_types=write_files)

    if write_files:
        list2file([rtp_text], os.path.join(ffpath, f"{hybrid_name}.rtp"))
        list2file([mtp_text], f"./{hybrid_name}.mtp")
        if bCharmm:
            list2file([_ARN_CONTENT], os.path.join(ffpath, f"{hybrid_name}.arn"))
        fn_restypes = pmx_data_file("mutff/residuetypes.dat")
        existing = _get_residue_types(fn_restypes)
        mode = "a" if os.path.isfile(fn_restypes) else "w"
        with open(fn_restypes, mode) as f:
            if hybrid_name not in existing:
                f.write(f"{hybrid_name:<6} Protein\n")

    if verbose:
        print(f"log_> covalent.mutate: hybrid {hybrid_name} "
              f"({residue1.resnA} -> {target})")
    return hybrid_name
