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

import collections
import glob
import logging
import os
import re
import warnings
from copy import deepcopy
from typing import Optional, Tuple

import numpy as np

from .. import library
from ..model import Model
from ..ffparser import RTPParser, load_ff_rtp
from ..utils import get_ff_path, list2file
from ..library import ff_data_file, pmx_data_file

from .mapping import (map_residues, attach_bonds_from_rtp, _element_of,
                      _canon_name, _heavy_neighbors, _is_h)
from .target import TargetResidue, complete_target_from_structure

# Reuse the battle-tested emit + ITP-fragment helpers from nsaa (import, do not
# modify).
from ..nsaa import (
    _merge_molecules, _make_bstate_dummies, _make_transition_dics, _unique_dummy_name,
    _update_bond_lists, _improps_as_atoms, _generate_dihedral_entries,
    _generate_improp_entries, _make_rtp, _make_mtp, _write_atp_fnb,
    _rtp_terminal, _assign_rtp_entries, _assign_mass_atp,
    _validate_atom_parameters, _get_residue_types, _ARN_CONTENT,
    _parse_itp_atoms, _parse_itp_bonds, _find_methyl_cap,
    _prepare_nsaa_fragment, _assign_itp_entries,
)


# ── coordinate placement ──────────────────────────────────────────────────────


logger = logging.getLogger(__name__)

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


def _h_count(res, i) -> int:
    return sum(1 for b in getattr(res.atoms[i], "bonds", [])
               for o in res.atoms if o is b and _is_h(o))


def _is_stereocentre(res, i) -> bool:
    """Four substituents, all distinguishable.

    Substituent COUNT alone is not enough.  Requiring only three mapped heavy
    neighbours fired on every planar sp2 carbon -- ASN ``CG``, ASP ``CG``, GLN
    ``CD`` are amide and carboxyl carbons with three substituents and no fourth --
    and on centres whose branches are equivalent, such as VAL ``CB`` with its two
    identical methyls.  Distinctness is judged by each substituent's
    (element, heavy degree, attached hydrogens) signature, which separates ILE's
    ``CG1`` (C,2,2) from ``CG2`` (C,1,3) while collapsing VAL's two methyls.
    """
    nb = _heavy_neighbors(res, i)
    nh = _h_count(res, i)
    if len(nb) + nh != 4 or len(nb) < 3:
        return False
    sig = [(_element_of(res.atoms[p]), len(_heavy_neighbors(res, p)),
            _h_count(res, p)) for p in nb]
    sig += [("H", 0, 0)] * nh
    return len(set(sig)) == 4


_CHIRAL_REF = ("N", "C", "CB")      # with CA as the centre


def _signed_volume(res, centre, three):
    """Triple product of *centre*'s bonds to *three*, or None if an atom is absent."""
    ix = {_canon_name(a.name): i for i, a in enumerate(res.atoms)}
    if centre not in ix or any(n not in ix for n in three):
        return None
    c = np.array(list(res.atoms[ix[centre]].x), float)
    v = [np.array(list(res.atoms[ix[n]].x), float) - c for n in three]
    return float(np.dot(np.cross(v[0], v[1]), v[2]))


def _match_chirality(residue1, residue2, pairs) -> bool:
    """Reflect the target when its CA configuration mirrors the WT's.

    An RTP is a bond graph and carries no stereochemistry, so the embedding in
    ``TargetResidue.from_rtp`` picks a configuration blind to L versus D: measured
    over a nine-residue peptide, SIX of nine came out as the D-enantiomer.  Only
    APPEARING atoms take their coordinates from the target, so a mirrored template
    is merely a poor starting geometry for a plain branch -- but for a stereocentre
    built out of appearing atoms it is the wrong diastereomer, and minimisation
    cannot repair that (the two are separated by a bond-breaking barrier).

    The WT is the reference because it comes from a real structure and is therefore
    L by construction.  Comparing against it needs no convention table and holds for
    any force field or atom-naming scheme.

    Returns True if the target was reflected.
    """
    s1 = _signed_volume(residue1, "CA", _CHIRAL_REF)
    s2 = _signed_volume(residue2, "CA", _CHIRAL_REF)
    if s1 is None or s2 is None:
        # glycine on either side: no CB, so CA is not a centre and there is no
        # reference triple.  A GLY WT leaves the target unchecked -- see the
        # module docstring note.
        return False
    if abs(s1) < 1e-6 or abs(s2) < 1e-6 or (s1 > 0) == (s2 > 0):
        return False

    # Reflect in a plane through the target's centroid: bond lengths and angles
    # are preserved, only the configuration flips.
    pts = np.array([list(a.x) for a in residue2.atoms], float)
    pts[:, 2] = 2.0 * float(pts[:, 2].mean()) - pts[:, 2]
    for a, row in zip(residue2.atoms, pts):
        a.x = [float(row[0]), float(row[1]), float(row[2])]
    logger.info("mutate: reflected %s -- its built CA configuration was the mirror "
                "of %s's (signed volume %+.2f vs %+.2f)",
                getattr(residue2, "resname", "?"),
                getattr(residue1, "resname", "?"), s2, s1)

    # A reflection inverts EVERY centre, so any other one that corresponds to a WT
    # centre may now disagree.  Report it rather than silently trading one wrong
    # configuration for another.
    mp = dict(pairs)
    for i, j in pairs:
        if _is_h(residue1.atoms[i]) or _canon_name(residue1.atoms[i].name) == "CA":
            continue
        if not _is_stereocentre(residue1, i):
            continue
        nb1 = [p for p in _heavy_neighbors(residue1, i) if p in mp]
        if len(nb1) < 3:
            continue
        n1 = [_canon_name(residue1.atoms[p].name) for p in nb1[:3]]
        n2 = [_canon_name(residue2.atoms[mp[p]].name) for p in nb1[:3]]
        a = _signed_volume(residue1, _canon_name(residue1.atoms[i].name), n1)
        b = _signed_volume(residue2, _canon_name(residue2.atoms[j].name), n2)
        if a is None or b is None or abs(a) < 1e-6 or abs(b) < 1e-6:
            continue
        if (a > 0) != (b > 0):
            logger.warning(
                "mutate: %s %s is now the mirror of %s %s after the CA correction; "
                "this residue has a second stereocentre and the two cannot both be "
                "satisfied by reflecting. Supply the target geometry explicitly "
                "(structure=/itp=) if the configuration matters.",
                getattr(residue2, "resname", "?"), residue2.atoms[j].name,
                getattr(residue1, "resname", "?"), residue1.atoms[i].name)
    return True


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
    # base RTP + any topic RTPs (nos.rtp, caps.rtp, ...); see ffparser.load_ff_rtp
    return load_ff_rtp(ff), get_ff_path(ff)


def _inter_residue_bonds(rtp, resname, term, ff=None):
    entr = (_rtp_terminal(rtp[resname], term, ff=ff, resname=resname)
            if term else rtp[resname])
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
                         term_b: Optional[str] = "__same__",
                         hybrid_name: Optional[str] = None,
                         write_ff_types: bool = True,
                         reverse: bool = False) -> Tuple[str, str, str]:
    """Build a hybrid residue from two parameterized residues via graph mapping.

    Both residues must already carry atom types/charges/masses and bond
    connectivity (``TargetResidue.from_rtp`` provides this for the target; a WT
    residue can be prepared with :func:`_assign_rtp_entries` +
    :func:`attach_bonds_from_rtp`).

    *term* patches state A's RTP bonded terms; *term_b* patches state B's (it
    defaults to *term*, so a plain mutation passes only *term*).  A terminal
    deletion's neighbour is the asymmetric case ``term=None, term_b='C'|'N'`` --
    an internal residue in state A that becomes a terminus in state B.
    *hybrid_name* overrides the generated ``X2Y`` name (e.g. ``VdeC``).

    Returns ``(rtp_text, mtp_text, hybrid_name)``.  If *write_ff_types*, new
    ``DUM_`` atom types are appended to the force field's ``atomtypes.atp`` /
    ``ffnonbonded.itp``.
    """
    if term_b == "__same__":
        term_b = term
    rtp, ffpath = _load_rtp(ff)
    bCharmm = "charmm" in ff

    _validate_atom_parameters(residue1, "state A")
    _validate_atom_parameters(residue2, "state B")
    _assign_mass_atp(residue1, residue2,
                     os.path.join(ffpath, "atomtypes.atp"), itp_file=None)

    # --- graph-MCS mapping (replaces domapping + polar-H patch) ---
    mapping = map_residues(residue1, residue2)
    pairs1, pairs2 = mapping.pairs1(), mapping.pairs2()

    # force the target's CA configuration to match the WT's before placing it
    _match_chirality(residue1, residue2, list(zip(pairs1, pairs2)))

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

    entr1 = _rtp_terminal(rtp[residue1.resname], term, ff=ff,
                          resname=residue1.resname) if term \
        else rtp[residue1.resname]
    if residue2.resname in rtp:
        entr2 = _rtp_terminal(rtp[residue2.resname], term_b, ff=ff,
                              resname=residue2.resname) if term_b \
            else rtp[residue2.resname]
    else:
        entr2 = {"atoms": [], "bonds": [], "diheds": [], "improps": [], "cmap": []}

    dih1 = _improps_as_atoms(entr1["diheds"], residue1, residue2)
    dih2 = _improps_as_atoms(entr2["diheds"], residue1, residue2, use_b=True)
    im1 = _improps_as_atoms(entr1["improps"], residue1, residue2)
    im2 = _improps_as_atoms(entr2["improps"], residue1, residue2, use_b=True)

    # reverse the alchemical direction: the target becomes state A and the WT
    # state B, so the emitted hybrid morphs product -> reactant.  Swapping the
    # A/B improper & dihedral lists flips which state each interaction belongs to
    # (the per-atom A/B fields are swapped just below, before emit).
    if reverse:
        dih1, dih2 = dih2, dih1
        im1, im2 = im2, im1

    cmap = rtp[residue1.resname].get("cmap", []) if bCharmm else []
    dihi_list = _generate_dihedral_entries(dih1, dih2, residue1, atom_pairs)
    ii_list = _generate_improp_entries(im1, im2, residue1)

    bond_neigh = _inter_residue_bonds(rtp, residue1.resnA if hasattr(
        residue1, "resnA") else residue1.resname, term, ff=ff)

    if reverse:
        # swap each atom's A and B state (types, charges, masses, names).  After
        # this residue1 carries the product as state A; _make_rtp/_make_mtp and
        # the model.write() of the mutant structure both read these swapped
        # fields, so the emitted hybrid and the output coordinates are consistent.
        for a in residue1.atoms:
            a.atomtype, a.atomtypeB = a.atomtypeB, a.atomtype
            a.q, a.qB = a.qB, a.q
            a.m, a.mB = a.mB, a.m
            a.name, a.nameB = a.nameB, a.name
            a.type, a.typeB = a.typeB, a.type
        # A WT-only atom (removed in the forward morph) carries the sentinel
        # name '<name>.gone' as its B-state name; after the swap that lands in
        # state A, which is not a valid structure-atom name.  Give it a clean
        # D-prefixed dummy name (pmx's added-atom convention), keeping RTP and
        # written coordinates consistent so pdb2gmx sorts it.  The name must stay
        # unique within the residue: 'D' + a 4-character name truncated to 4
        # would map ARG's HH11/HH12 both onto DHH1 (and LEU's HD11/HD12/HD13
        # onto DHD1), and pdb2gmx then keeps only one of the pair, silently
        # dropping atoms and their charge from state B.
        # Hydrogens get pmx's HV<n> counter, exactly as the FORWARD generator does
        # (nsaa._merge_molecules) and exactly as the shipped library spells them: the
        # pre-built ALA->ARG hybrid A2R names its nine hydrogen dummies HV1..HV9.
        # Using the same convention here makes a generated reverse hybrid come out
        # atom-for-atom identical to the library entry for the same morph, which is
        # the point -- reversing X->Y must be indistinguishable from asking for Y->X.
        # It also removes the collision class outright: every name that used to clash
        # under ('D'+name)[:4] was a hydrogen (ARG HH11/HH12, LEU HD11/HD12/HD13,
        # VAL/ILE/THR HG1x, ASN HD2x, GLN HE2x), and an HV counter cannot collide.
        # Heavy atoms keep 'D'+name, still guarded for the 4-character cap.
        used = {a.name for a in residue1.atoms if not a.name.endswith('.gone')}
        hv = 1
        for a in residue1.atoms:
            if not a.name.endswith('.gone'):
                continue
            base = a.name[:-len('.gone')]
            if base.startswith('H') or (base[0].isdigit() and len(base) > 1
                                        and base[1] == 'H'):
                while ('HV%d' % hv) in used:
                    hv += 1
                a.name = 'HV%d' % hv
                hv += 1
            else:
                a.name = _unique_dummy_name(base, used)
            used.add(a.name)

    _auto_named = hybrid_name is None
    if hybrid_name is None:
        # Central namer (library.hybrid_name).  Internal: L2<first+last of target>
        # so prefix-sharing PTM targets stay distinct (MLZ->K2MZ vs MLY->K2MY vs
        # M3L->K2ML).  Terminal (term set): the classic library spelling plus a
        # lowercase terminus letter (L2Yn / L2Tc) -- a terminal hybrid is a
        # different residue from the internal one of the same morph, so it gets
        # its own name instead of shadowing the library's L2Y.  A reverse hybrid
        # front-loads the target so it never collides with its forward
        # counterpart in the shared nsaa.rtp/nsaa.mtp.
        hybrid_name = library.hybrid_name(residue1.resname, mut_resname,
                                          reverse=reverse, term=term)
    if reverse:
        residue1.resnA, residue1.resnB = mut_resname, residue1.resname
    else:
        residue1.resnA, residue1.resnB = residue1.resname, mut_resname
    if _auto_named and write_ff_types and _library_owns(ffpath, hybrid_name)[0]:
        # The generated name is one the curated library already uses.  That is
        # only a problem when OUR hybrid differs from theirs -- an identical one
        # may keep the spelling and reuse the shipped entry.  _make_mtp is pure,
        # so probe with it before committing to a name.
        _probe = "".join(_make_mtp(residue1, ii_list, dihi_list, rotations=None))
        if _differs_from_curated(ffpath, hybrid_name, _probe):
            hybrid_name = _generative_name(ffpath, hybrid_name,
                                           (residue1.resnA, residue1.resnB))
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
    (geometry override, parameters from the RTP) > RTP name (default: parameters +
    RDKit-embedded geometry from the force-field RTP entry).
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
        # If the structure is MISSING atoms that the RTP defines (e.g. a
        # heavy-atom-only crystal template with no hydrogens), complete the atom
        # set from the RTP and overlay only the supplied heavy-atom coordinates.
        # Complete structures are left untouched, so every other structure=
        # workflow is unchanged.
        residue2 = complete_target_from_structure(residue2, target, ff, term)
        _assign_rtp_entries(residue2, rtp, term, ff=ff)
        if not any(getattr(a, "bonds", None) for a in residue2.atoms):
            attach_bonds_from_rtp(residue2, ff)
        return residue2

    # default: parameters + embedded geometry from the RTP entry.  When the
    # residue sits at a terminus, build the target in its terminal-patched form
    # too -- otherwise the WT (terminal) morphs to an internal target and the
    # terminus charge (e.g. the N-terminal +1) is lost in state B.
    return TargetResidue.from_rtp(target, ff, term=term).residue


# ── hybrid accumulation (single per-force-field nsaa.rtp / nsaa.mtp) ───────────

def _entry_morph(text, name):
    """Return the (resnA, resnB) recorded in a ``[ name ] ; A -> B`` header.

    Both _make_rtp/_make_mtp stamp each entry with its source->target morph, so
    this lets _accumulate_hybrid tell a rebuild of the *same* hybrid (replace)
    from a *different* morph that hashed to the same 4-char name (collision).
    Returns None when the name is absent or carries no morph comment.
    """
    m = re.search(r'^\[\s*' + re.escape(name) + r'\s*\]\s*;\s*(\S+)\s*->\s*(\S+)',
                  text, re.MULTILINE)
    return (m.group(1), m.group(2)) if m else None


def _library_owns(ffpath, name):
    """-> ((resnA, resnB), basename) if a CURATED rtp (not nsaa.*) defines *name*."""
    for other in sorted(glob.glob(os.path.join(ffpath, "*.rtp"))):
        if os.path.basename(other).startswith("nsaa."):
            continue
        morph = _entry_morph(open(other, errors="replace").read(), name)
        if morph:
            return morph, os.path.basename(other)
    return None, None


def _split_named_blocks(text):
    """-> {NAME: block_text} for a multi-entry mtp/rtp.

    Top-level entries are spelled ``[ NAME ]`` in upper case/digits; the
    subsections inside one (``[ atoms ]``, ``[ bonds ]``, ...) are lower case, so
    an upper-case-only lookahead splits entries without splitting their guts.
    """
    out = {}
    for part in re.split(r'\n(?=\[\s*[A-Z0-9]+\s*\])', text):
        m = re.match(r'\[\s*([A-Z0-9]+)\s*\]', part.strip())
        if m:
            out[m.group(1)] = part
    return out


def _mtp_atom_rows(block):
    """(typeA, qA, mA, typeB, qB, mB) per atom of an mtp ``[ atoms ]`` section.

    Atom NAMES are excluded on purpose: the curation and the engine may pick
    different members of a set of equivalent hydrogens and spell H names
    differently (1HB vs HB1), which is not a chemical difference.
    """
    rows, inseg = [], False
    for line in block.splitlines():
        t = line.split(';')[0].strip()
        if t.startswith('[ atoms'):
            inseg = True
            continue
        if inseg and t.startswith('['):
            break
        f = t.split()
        if inseg and len(f) >= 8:
            try:
                rows.append((f[1], round(float(f[2]), 4), round(float(f[4]), 3),
                             f[5], round(float(f[6]), 4), round(float(f[7]), 3)))
            except ValueError:
                pass
    return collections.Counter(rows)


def _differs_from_curated(ffpath, name, mtp_text):
    """True when the hybrid just built is NOT the curated entry of the same name.

    Only then does it need a name of its own: an identical hybrid can safely keep
    the library's spelling and reuse the shipped entry, which is what
    _accumulate_hybrid already does.
    """
    mine = _mtp_atom_rows(mtp_text)
    for other in sorted(glob.glob(os.path.join(ffpath, "*.mtp"))):
        if os.path.basename(other).startswith("nsaa."):
            continue
        block = _split_named_blocks(open(other, errors="replace").read()).get(name)
        if block is not None:
            return _mtp_atom_rows(block) != mine
    # the name is curated in an rtp but has no mtp entry to compare against --
    # cannot prove the two are the same, so do not risk shadowing it
    return True


def _generative_name(ffpath, name, morph):
    """A name of its own for a generative hybrid that differs from the library.

    The forward namer deliberately reproduces the library's spelling, so a
    generative build of a *canonical* morph lands on the name mutres.rtp already
    uses.  Where the two are genuinely different objects -- a different mapping,
    a different dummy count -- pdb2gmx would glob ``*.rtp``, pick the curated
    entry and leave the written structure unmatched by it (``Atom HV1 in residue
    C2S was not found in rtp entry C2S``).  Mark ours with a trailing ``G``
    instead, inside the 4-character residue-name cap: ``C2S`` -> ``C2SG``.

    A candidate already holding THIS morph in nsaa.rtp is reused, so rebuilding
    the same hybrid is idempotent; one held by a DIFFERENT morph is skipped
    (``CM2I`` cannot take ``CM2G``, which is the library's CYM->GLY).
    """
    up = lambda t: tuple(x.upper() for x in t)
    cands = (["%sG" % name] if len(name) < 4 else []) + \
            [name[:3] + "G", name[:2] + "G" + name[-1], "G" + name[:3]]
    nsaa = os.path.join(ffpath, "nsaa.rtp")
    have = open(nsaa, errors="replace").read() if os.path.isfile(nsaa) else ""
    for c in cands:
        if _library_owns(ffpath, c)[0]:
            continue
        mine = _entry_morph(have, c)
        if mine is None or up(mine) == up(morph):
            logger.info("generative %s->%s differs from the curated %s; naming it "
                        "%s so it does not shadow the library entry",
                        morph[0], morph[1], name, c)
            return c
    raise ValueError(
        "cannot find a free generative name for %s->%s near '%s' (tried %s); "
        "pass an explicit hybrid_name=" % (morph[0], morph[1], name,
                                           ", ".join(cands)))


def _accumulate_hybrid(ffpath, hybrid_name, rtp_text, mtp_text, bCharmm):
    """Upsert an engine-built hybrid into the force field's nsaa.{rtp,mtp}.

    All generative hybrids -- point mutations, PTMs, cross-link endpoints and
    terminal-deletion neighbours alike -- accumulate into single per-force-field
    files (mirroring the mutres.rtp/mutres.mtp convention) rather than scattering
    one file per hybrid.  pdb2gmx globs ``*.rtp`` in the FF dir and gen_hybrid_top
    reads multi-entry MTPs, so a single file serves both.  Re-building a hybrid
    replaces its entry.
    """
    # Guard the 4-char residue-name cap: refuse to overwrite an existing entry
    # that belongs to a *different* morph (silent overwrite would corrupt a
    # system holding both, e.g. a mono- and a di-methyl-Lys).  Rebuilding the
    # same morph is fine (it replaces its own entry).
    rtp_path = os.path.join(ffpath, "nsaa.rtp")
    incoming = _entry_morph(rtp_text, hybrid_name)
    # A reverse hybrid of a canonical morph is now named exactly as the shipped
    # library names that mutation (reverse of ARG->ALA is 'A2R'), so the entry we are
    # about to write may ALREADY exist in mutres.rtp.  Writing a second definition
    # under the same name would leave two [ A2R ] blocks in one force field -- pdb2gmx
    # globs *.rtp -- which is the nsaa-pollution failure mode.  The shipped entry is
    # hand-curated, so it wins: skip the write.  A same-named entry for a DIFFERENT
    # morph is a genuine collision and still raises.
    if incoming:
        for other in sorted(glob.glob(os.path.join(ffpath, "*.rtp"))):
            if os.path.basename(other).startswith("nsaa."):
                continue
            shipped = _entry_morph(open(other, errors="replace").read(), hybrid_name)
            if not shipped:
                continue
            if tuple(x.upper() for x in shipped) == tuple(x.upper() for x in incoming):
                logger.info("hybrid %s already shipped in %s (%s->%s); keeping the "
                            "library entry", hybrid_name, os.path.basename(other),
                            *shipped)
                return
            raise ValueError(
                f"hybrid name '{hybrid_name}' is already used by "
                f"{os.path.basename(other)} for {shipped[0]}->{shipped[1]}; refusing "
                f"to shadow it with {incoming[0]}->{incoming[1]}. Pass an explicit "
                f"hybrid_name= to disambiguate.")
    if incoming and os.path.isfile(rtp_path):
        existing = _entry_morph(open(rtp_path).read(), hybrid_name)
        if existing and existing != incoming:
            raise ValueError(
                f"hybrid name '{hybrid_name}' already maps "
                f"{existing[0]}->{existing[1]} in {rtp_path}; refusing to reuse "
                f"it for {incoming[0]}->{incoming[1]}. Pass an explicit "
                f"hybrid_name= to disambiguate.")
    rtp_header, rtp_entry = _split_rtp_header(rtp_text, hybrid_name)
    _upsert_entry(os.path.join(ffpath, "nsaa.rtp"), hybrid_name, rtp_entry,
                  header=rtp_header)
    _upsert_entry(os.path.join(ffpath, "nsaa.mtp"), hybrid_name, mtp_text)
    if bCharmm:
        # generic atom-renaming spec (constant content) -- write once
        arn = os.path.join(ffpath, "nsaa.arn")
        if not os.path.isfile(arn):
            list2file([_ARN_CONTENT], arn)
    fn_restypes = ff_data_file(ffpath, "residuetypes.dat")
    existing = _get_residue_types(fn_restypes)
    if hybrid_name not in existing:
        # residuetypes.dat may lack a trailing newline; GROMACS rejects a
        # merged line, so ensure the previous entry is terminated first.
        need_nl = (os.path.isfile(fn_restypes)
                   and os.path.getsize(fn_restypes) > 0
                   and open(fn_restypes, "rb").read()[-1:] != b"\n")
        with open(fn_restypes, "a" if os.path.isfile(fn_restypes) else "w") as f:
            f.write(("\n" if need_nl else "") + f"{hybrid_name:<6} Protein\n")


# ── model-level mutate ────────────────────────────────────────────────────────

def _apply_ref_geometry(residue2, ref_path, mut_resid, mut_chain):
    """Override the target residue's coordinates from a reference (state-B)
    structure's residue at the same position, matched by atom name.

    This is the generative-path equivalent of the classic ``-fB``: the modified
    residue's B-state geometry follows a supplied structure of the modified
    state rather than the RDKit-embedded / fragment geometry.  Atoms absent from
    the reference keep their built geometry.
    """
    mB = Model(ref_path, bPDBTER=True, rename_atoms=True, scale_coords='A')
    resB = mB.fetch_residue(idx=mut_resid, chain=mut_chain)
    byname = {a.name: a for a in resB.atoms}
    for a in residue2.atoms:
        if a.name in byname:
            a.x = list(byname[a.name].x)



def mutate(model, mut_resid: int, target: str, ff: str,
           structure: Optional[str] = None, itp: Optional[str] = None,
           refB: Optional[str] = None,
           mut_chain: Optional[str] = None, inplace: bool = False,
           write_files: bool = True, verbose: bool = False,
           reverse: bool = False, hybrid_name: Optional[str] = None):
    """Mutate ``mut_resid`` in *model* to *target*, via one of three front ends.

    Parameters
    ----------
    target : str
        Target residue name (also used for the hybrid name and RTP lookup).
    hybrid_name : str, optional
        Override the generated 4-character hybrid name.  ``library.hybrid_name``
        builds it from the source one-letter code plus the FIRST and LAST letters
        of the target's 3-letter code, which is **not** collision-free: THR and
        TYR both yield ``L2TR`` from LEU.  The cross-morph guard in
        :func:`_accumulate_hybrid` refuses to reuse a name for a different morph
        and tells the caller to pass this -- so it has to be reachable from here.
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
    # Accept the classic one-letter residue shorthand for the target and resolve
    # it to the 3-letter RTP name the generative engine needs.  The classic
    # library path resolves these itself (via _convert_aa_name); the generative
    # path did not, so a one-letter target at a terminus / PTM / NSAA -- e.g.
    # `pmx mutate --script '1 V'` on an N-terminal residue -- failed with
    # "Residue 'V' not found in the RTP".
    if len(target) == 1:
        target = library._aacids_dic.get(target.upper(), target)

    rtp, ffpath = _load_rtp(ff)
    bCharmm = "charmm" in ff

    m = model if inplace else deepcopy(model)
    if m.unity == "nm":
        m.nm2a()

    residue1 = m.fetch_residue(idx=mut_resid, chain=mut_chain)
    chain = residue1.chain
    if residue1 == chain.cterminus() and residue1 == chain.nterminus():
        _reject_lone_residue(residue1, "covalent.mutate")
    term = ("C" if residue1 == chain.cterminus()
            else ("N" if residue1 == chain.nterminus() else None))

    # parameterize + bond the WT residue from the force field
    _assign_rtp_entries(residue1, rtp, term, ff=ff)
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

    # optional state-B reference (classic -fB): take the modified residue's
    # geometry from a supplied structure of the mutated state at this position.
    if refB is not None:
        _apply_ref_geometry(residue2, refB, mut_resid, mut_chain)

    rtp_text, mtp_text, hybrid_name = build_hybrid_residue(
        residue1, residue2, ff, target, term=term, write_ff_types=write_files,
        reverse=reverse, hybrid_name=hybrid_name)

    if write_files:
        _accumulate_hybrid(ffpath, hybrid_name, rtp_text, mtp_text, bCharmm)

    if verbose:
        print(f"log_> covalent.mutate: hybrid {hybrid_name} "
              f"({residue1.resnA} -> {residue1.resnB}"
              f"{', reverse' if reverse else ''}); "
              f"MTP -> {os.path.join(ffpath, 'nsaa.mtp')}")
    return hybrid_name


# ── terminal deletion (a residue morph through the same core) ──────────────────

def _reject_lone_residue(residue, op: str) -> None:
    """Refuse an op on a chain of ONE residue (it is N- AND C-terminal).

    Such a residue needs BOTH terminal patches, and the engine applies only one
    (``term`` resolves to a single side), so it would silently emit a hybrid whose
    backbone is patched on one end and internal on the other.  A **dipeptide is
    the minimum supported input**; a lone residue is out of scope and is rejected
    here rather than half-handled.  The matching case on the deletion side -- a
    deletion whose survivor would be a lone residue -- is caught by the strict
    ``_validate_atom_parameters`` call in :func:`delete_terminal`.
    """
    chain = residue.chain
    if chain is None:
        return
    if len(chain.residues) < 2:
        raise ValueError(
            f"{op}: residue {residue.resname}{residue.id} is the only residue in "
            f"chain '{getattr(chain, 'id', '?')}', so it is both the N- and the "
            f"C-terminus and would need both terminal patches.  pmx applies one, "
            f"so this is refused: a dipeptide is the minimum supported input.")


def _term_del_code(resname: str) -> str:
    """One-letter code for the deletion-neighbour hybrid name (``V`` -> ``VdeC``).

    Uses ``library._ext_one_letter`` -- the canonical map that already
    disambiguates the isosteric residue-name collisions (HSD->H, HSE->X,
    HSP->Z, CYM->CM), and from which ``library._term_deletion_residues`` (the
    recognised XdeC/XdeN protein-residue set) is derived -- so a generated name
    is guaranteed to be recognised and never collides across variants.
    """
    return library._ext_one_letter.get(resname,
                                        library._one_letter.get(resname,
                                                                resname[0]))


def delete_terminal(model, del_resid: int, ff: str,
                    del_chain: Optional[str] = None, inplace: bool = False,
                    write_files: bool = True, verbose: bool = False):
    """Alchemically delete a terminal residue, composed from the core engine.

    A terminal deletion is two residue morphs plus a boundary decouple:

    * the **neighbour** (n-1 for a C-terminal deletion, n+1 for an N-terminal)
      morphs from an internal residue to a terminus -- an ordinary
      ``build_hybrid_residue`` graph map of ``internal -> terminal`` (the
      terminal state built by ``TargetResidue.from_rtp(..., term=...)``), so it
      grows the terminus cap atom (``OT2`` / extra amide H's) as an appearing
      dummy and retypes the backbone.  Emitted as the ``XdeC``/``XdeN`` hybrid.
    * the **target** residue becomes all-dummy in state B (done at gen_hybrid_top,
      which also decouples the boundary dihedrals).

    This replaces the bespoke ``alchemy.apply_terminal_deletion`` +
    ``scripts.generate_term_deletion``: the neighbour hybrid is now generated on
    the fly into ``nsaa.{rtp,mtp}`` like every other engine hybrid, rather than
    pre-built per force field.
    """
    rtp, ffpath = _load_rtp(ff)
    bCharmm = "charmm" in ff

    m = model if inplace else deepcopy(model)
    if m.unity == "nm":
        m.nm2a()

    del_res = m.fetch_residue(idx=del_resid, chain=del_chain)
    chain = del_res.chain
    residues = chain.residues
    idx = residues.index(del_res)

    is_c = del_res == chain.cterminus()
    is_n = del_res == chain.nterminus()
    if not (is_c or is_n):
        raise ValueError(
            "Terminal deletion requested but residue %s (%s) is not at a "
            "terminus" % (del_res.id, del_res.resname))
    if is_c:
        if idx == 0:
            raise ValueError("Cannot delete C-terminal residue: it is the only "
                             "residue in the chain")
        neighbor, term, suffix = residues[idx - 1], "C", "deC"
    else:
        if idx >= len(residues) - 1:
            raise ValueError("Cannot delete N-terminal residue: it is the only "
                             "residue in the chain")
        neighbor, term, suffix = residues[idx + 1], "N", "deN"

    # parameterize + bond the neighbour as an internal residue (state A)
    _assign_rtp_entries(neighbor, rtp, None)
    # The neighbour is typed as an INTERNAL residue in state A.  If it is also the
    # opposite terminus (a one-residue chain survives the deletion) it still
    # carries that terminus's atoms -- OT1/OT2 or H1/H2/H3 -- which the internal
    # entry does not define, so they would stay untyped and surface as an
    # unparseable three-field nsaa.rtp row inside pdb2gmx.  Fail here instead.
    _validate_atom_parameters(neighbor, f"neighbour of a {term}-terminal deletion",
                              strict=True)
    if not any(getattr(a, "bonds", None) for a in neighbor.atoms):
        attach_bonds_from_rtp(neighbor, ff)

    # state B: the same residue in its terminal-patched form (embedded geometry,
    # superposed onto the neighbour's real coords -- the cap dummy rides along)
    residue2 = TargetResidue.from_rtp(neighbor.resname, ff, term=term).residue

    hybrid_name = _term_del_code(neighbor.resname) + suffix
    rtp_text, mtp_text, hn = build_hybrid_residue(
        neighbor, residue2, ff, neighbor.resname, term=None, term_b=term,
        hybrid_name=hybrid_name, write_ff_types=write_files)

    # The deleted residue is a WT -> all-dummy morph: every atom disappears in
    # state B.  gen_hybrid_top applies that dummification structurally, but the
    # backbone DUM_ types it introduces (e.g. DUM_NH1, never needed by an
    # ordinary sidechain mutation) must be registered in the force field.  Build
    # the all-dummy hybrid through the engine so _make_bstate_dummies +
    # _write_atp_fnb register exactly those types -- the same mechanism every
    # engine hybrid uses.
    if write_files:
        # ff= so the terminal patch comes from this force field's .tdb: the
        # DUM_ types registered here must be the ones pdb2gmx will actually
        # assign to the deleted residue (amber N-term CA is CX, not CHARMM CT1).
        _assign_rtp_entries(del_res, rtp, term, ff=ff)
        # strict: an untyped atom here means the survivor needs a patch this
        # entry does not provide (e.g. a one-residue chain after the deletion is
        # both N- and C-terminal).  Fail now, named, not later inside pdb2gmx.
        _validate_atom_parameters(del_res, f"{term}-terminal deletion target",
                                  strict=True)
        _make_bstate_dummies(del_res)
        _write_atp_fnb("atomtypes.atp", "ffnonbonded.itp", del_res, ffpath)

    if write_files:
        _accumulate_hybrid(ffpath, hn, rtp_text, mtp_text, bCharmm)

    if verbose:
        print(f"log_> covalent.delete_terminal: {term}-terminal deletion of "
              f"residue {del_res.id} ({del_res.resname}); neighbour "
              f"{neighbor.id} -> hybrid {hn}")
    return m if not inplace else hn


def amidate_cterminus(model, resid: int, ff: str, chain: Optional[str] = None,
                      inplace: bool = False, write_files: bool = True,
                      verbose: bool = False, reverse: bool = False,
                      hybrid_name: Optional[str] = None):
    """Alchemically amidate a C-terminus: ``-COO-`` (A) -> ``-C(=O)NH2`` (B).

    The experimental "amidation" that blocks the free C-terminal carboxylate (a
    C-degron recognition test).  It is charge-changing (-1 -> 0): the second
    carboxylate oxygen disappears and an amide NH2 appears -- built through the
    same hybrid-residue engine as a mutation, via the CHARMM CT2 terminal patch
    (``term_b='C_amide'``).  A *terminal modification*, the addition companion to
    :func:`delete_terminal`.
    """
    rtp, ffpath = _load_rtp(ff)
    bCharmm = "charmm" in ff
    m = model if inplace else deepcopy(model)
    if m.unity == "nm":
        m.nm2a()

    res = m.fetch_residue(idx=resid, chain=chain)
    _reject_lone_residue(res, "amidate_cterminus")
    if res != res.chain.cterminus():
        raise ValueError("amidate_cterminus: residue %s (%s) is not the "
                         "C-terminus" % (res.id, res.resname))
    # the residue already carries its carboxylate (OT1/OT2 from pdb2gmx), so it
    # must be typed with the C-terminal patch -- not the internal RTP, which
    # would leave OT2 untyped and empty out the hybrid.
    _assign_rtp_entries(res, rtp, "C", ff=ff)
    if not any(getattr(a, "bonds", None) for a in res.atoms):
        attach_bonds_from_rtp(res, ff)

    residue2 = TargetResidue.from_rtp(res.resname, ff, term="C_amide").residue
    code = library._ext_one_letter.get(res.resname, res.resname[:1])
    hybrid_name = hybrid_name or (code + "am")[:4]   # caller may override
    rtp_text, mtp_text, hn = build_hybrid_residue(
        res, residue2, ff, res.resname, term="C", term_b="C_amide", reverse=reverse,
        hybrid_name=hybrid_name, write_ff_types=write_files)
    if write_files:
        _accumulate_hybrid(ffpath, hn, rtp_text, mtp_text, bCharmm)
    if verbose:
        print(f"log_> covalent.amidate_cterminus: {res.id} ({res.resname}) "
              f"-COO- -> -CONH2  (hybrid {hn})")
    return m if not inplace else hn
