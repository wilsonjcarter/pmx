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
import warnings
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

# Reuse the battle-tested emit + ITP-fragment helpers from nsaa (import, do not
# modify).
from ..nsaa import (
    _merge_molecules, _make_bstate_dummies, _make_transition_dics,
    _update_bond_lists, _improps_as_atoms, _generate_dihedral_entries,
    _generate_improp_entries, _make_rtp, _make_mtp, _write_atp_fnb,
    _rtp_terminal, _assign_rtp_entries, _assign_mass_atp,
    _validate_atom_parameters, _get_residue_types, _ARN_CONTENT,
    _parse_itp_atoms, _parse_itp_bonds, _find_methyl_cap,
    _prepare_nsaa_fragment, _assign_itp_entries,
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


# section headers *within* a residue entry -- used to find where one entry ends
_SECTION_KW = {"atoms", "bonds", "angles", "impropers", "dihedrals", "cmap",
               "bondedtypes", "exclusions", "pairs", "morphes", "coords",
               "rotations", "settles", "virtual_sites"}


def _is_entry_header(line: str) -> bool:
    """True for a top-level residue header '[ NAME ]' (not a '[ atoms ]' etc.)."""
    s = line.strip()
    if not (s.startswith("[") and "]" in s):
        return False
    tok = s[1:s.index("]")].split()
    return bool(tok) and tok[0].lower() not in _SECTION_KW


def _split_rtp_header(rtp_text: str, name: str):
    """Split _make_rtp output into (file header, single residue entry)."""
    idx = rtp_text.find(f"[ {name} ]")
    return (rtp_text[:idx], rtp_text[idx:]) if idx != -1 else ("", rtp_text)


def _upsert_entry(path: str, name: str, entry_text: str, header: str = "") -> None:
    """Insert or replace the ``[ name ]`` block in a multi-entry RTP/MTP file.

    Mirrors the single-file mutres.rtp/mutres.mtp convention: one accumulating
    file per force field rather than a new file per hybrid.  Re-mutating the
    same hybrid replaces its entry instead of duplicating it.
    """
    if not entry_text.endswith("\n"):
        entry_text += "\n"
    if not os.path.isfile(path):
        with open(path, "w") as f:
            f.write(header + entry_text)
        return
    lines = open(path).read().splitlines(keepends=True)
    start = next((i for i, l in enumerate(lines)
                  if l.strip().startswith(f"[ {name} ]")), None)
    if start is None:
        with open(path, "a") as f:
            if lines and not lines[-1].endswith("\n"):
                f.write("\n")
            f.write("\n" + entry_text)
        return
    end = next((j for j in range(start + 1, len(lines))
                if _is_entry_header(lines[j])), len(lines))
    lines[start:end] = [entry_text]
    with open(path, "w") as f:
        f.writelines(lines)


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


# ── target-residue front ends ─────────────────────────────────────────────────

def _attach_itp_bonds(residue2, itp_file: str) -> None:
    """Give an ITP-parameterized residue connectivity (for the graph mapper)."""
    nr2name = {a["nr"]: a["name"] for a in _parse_itp_atoms(itp_file)}
    by_name = {a.name: a for a in residue2.atoms}
    for ai, aj in _parse_itp_bonds(itp_file):
        x, y = by_name.get(nr2name.get(ai)), by_name.get(nr2name.get(aj))
        if x is not None and y is not None and y not in x.bonds:
            x.bonds.append(y)
            y.bonds.append(x)


def _build_target_residue(target, ff, rtp, term, residue1,
                          structure=None, itp=None):
    """Build the (parameterized, bonded) target residue from the chosen front end.

    Precedence: ``itp`` (arbitrary parameterized fragment) > ``structure``
    (geometry override, parameters from the RTP) > RTP name (default:
    parameters + RDKit-embedded geometry from the force-field RTP entry).
    """
    if itp is not None:
        # NSAA capped-fragment workflow: parameters from the ITP; geometry from
        # the ITP bond graph (RDKit-embedded via TargetResidue.from_itp, exactly
        # like an RTP target) -- so no structure file is needed.  A structure may
        # still be passed to override the embedded geometry.  The backbone is
        # grafted from the WT (fragment is CH3-capped at CA).  Reuses pmx.nsaa's
        # fragment machinery.
        if structure is not None:
            model2 = Model()
            model2.read(structure)
            residue2 = model2.residues[0]
            residue2.set_resname(target)
        else:
            residue2 = TargetResidue.from_itp(itp, target).residue
        cap_C, cap_Hs, cap_CA, cap_CA_Hs = _find_methyl_cap(
            _parse_itp_atoms(itp), _parse_itp_bonds(itp))
        has_bb = all(any(a.name == n for a in residue2.atoms)
                     for n in ("N", "C", "O"))
        if not has_bb:
            if cap_C is None:
                raise RuntimeError(
                    "ITP fragment has no backbone N/C/O and no CH3 cap at CA.")
            _prepare_nsaa_fragment(residue2, residue1, cap_C, cap_Hs,
                                   ca_atom_name=cap_CA, ca_H_names=cap_CA_Hs,
                                   term=term)
        _assign_itp_entries(residue2, itp, rtp, backbone_resname=residue1.resname)
        _attach_itp_bonds(residue2, itp)
        return residue2

    if structure is not None:
        # geometry from the structure; parameters + connectivity from the RTP
        model2 = Model()
        model2.read(structure)
        residue2 = model2.residues[0]
        residue2.set_resname(target)
        _assign_rtp_entries(residue2, rtp, term)
        if not any(getattr(a, "bonds", None) for a in residue2.atoms):
            attach_bonds_from_rtp(residue2, ff)
        return residue2

    # default: parameters + embedded geometry from the RTP entry
    return TargetResidue.from_rtp(target, ff).residue


# ── model-level mutate ────────────────────────────────────────────────────────

def mutate(model, mut_resid: int, target: str, ff: str,
           structure: Optional[str] = None, itp: Optional[str] = None,
           mut_chain: Optional[str] = None, inplace: bool = False,
           write_files: bool = True, verbose: bool = False):
    """Mutate ``mut_resid`` in *model* to *target*, via one of three front ends.

    Parameters
    ----------
    target : str
        Target residue name (also used for the hybrid name and RTP lookup).
    structure : str, optional
        PDB/GRO with the target residue's geometry.  Overrides the default
        RDKit-embedded geometry; parameters still come from the RTP (or from
        *itp* if given).
    itp : str, optional
        GROMACS ITP fragment parameterizing the target (for residues not in the
        force-field RTP).  Requires *structure*.  Uses the NSAA capped-fragment
        workflow.

    With neither *structure* nor *itp*, the target is built entirely from the
    force-field RTP entry (parameters + embedded geometry).

    Accumulates the hybrid into the force field's ``nsaa.rtp``/``nsaa.mtp`` and
    returns the hybrid residue name.
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

    # Protonation check: a hydrogen-free WT residue (e.g. a raw crystal
    # structure) leaves the hybrid's backbone/polar hydrogens as dummies, so an
    # improper such as (N -C CA HN) references a now-absent atom and pdb2gmx
    # rejects the topology.  Warn up front rather than fail cryptically later.
    if not any(_element_of(a) == "H" for a in residue1.atoms):
        warnings.warn(
            f"WT residue {mut_resid} ({residue1.resname}) has no hydrogens -- "
            "the structure appears unprotonated. Its polar/backbone hydrogens "
            "will be treated as appearing dummies, which yields an invalid "
            "topology (pdb2gmx will reject a dangling improper). Protonate the "
            "structure first (e.g. run gmx pdb2gmx on the WT and mutate the "
            "protonated output).", stacklevel=2)

    # build the target residue from the selected front end (RTP / structure / itp)
    residue2 = _build_target_residue(target, ff, rtp, term, residue1,
                                     structure=structure, itp=itp)

    rtp_text, mtp_text, hybrid_name = build_hybrid_residue(
        residue1, residue2, ff, target, term=term, write_ff_types=write_files)

    if write_files:
        # Accumulate all engine-built hybrids into single per-force-field files
        # (nsaa.rtp / nsaa.mtp), mirroring the mutres.rtp/mutres.mtp convention,
        # instead of scattering one file per hybrid.  pdb2gmx globs *.rtp in the
        # FF dir, and gen_hybrid_top reads multi-entry MTPs -- so a single file
        # works for both.  Re-mutating a hybrid replaces its entry.
        rtp_header, rtp_entry = _split_rtp_header(rtp_text, hybrid_name)
        _upsert_entry(os.path.join(ffpath, "nsaa.rtp"), hybrid_name, rtp_entry,
                      header=rtp_header)
        _upsert_entry(os.path.join(ffpath, "nsaa.mtp"), hybrid_name, mtp_text)
        if bCharmm:
            # generic atom-renaming spec (constant content) -- write once
            arn = os.path.join(ffpath, "nsaa.arn")
            if not os.path.isfile(arn):
                list2file([_ARN_CONTENT], arn)
        fn_restypes = pmx_data_file("mutff/residuetypes.dat")
        existing = _get_residue_types(fn_restypes)
        if hybrid_name not in existing:
            # residuetypes.dat may lack a trailing newline; GROMACS rejects a
            # merged line, so ensure the previous entry is terminated first.
            need_nl = (os.path.isfile(fn_restypes)
                       and os.path.getsize(fn_restypes) > 0
                       and open(fn_restypes, "rb").read()[-1:] != b"\n")
            with open(fn_restypes, "a" if os.path.isfile(fn_restypes) else "w") as f:
                f.write(("\n" if need_nl else "") + f"{hybrid_name:<6} Protein\n")

    if verbose:
        print(f"log_> covalent.mutate: hybrid {hybrid_name} "
              f"({residue1.resnA} -> {target}); "
              f"MTP -> {os.path.join(ffpath, 'nsaa.mtp')}")
    return hybrid_name
