"""Auto-generated dual-topology (split-atom) hybrids for cross-links.

For a bond-forming covalent modification (disulfide `SG-SG'`, NOS `NZ-OD`, or any
cross-link), GROMACS generates non-bonded exclusions/1-4 pairs topologically, so
the forming bond would corrupt the free state.  The proven fix is to *decouple*
the reacting sidechain: duplicate it into two branches, one real in state A
(free), the other real in state B (bound), so everything the new bond touches
lives on the bound branch -- which is a dummy in state A, making its exclusions
harmless there.  This is the split-atom scheme of the shipped `C2CD` residue,
generated automatically instead of hand-built.

Step 1 (this module, so far): **walk the sidechain from the bond atom and
identify what to dummify** -- the reacting atoms to duplicate and the common
boundary atom where the two branches reconnect.
"""

from __future__ import annotations

import os
from typing import Dict, List, Optional, Tuple

from .. import library
from ..atom import Atom
from .mapping import _element_of

# backbone atoms are never duplicated -- they are the common scaffold
_BACKBONE = {"N", "HN", "H", "CA", "HA", "C", "O", "OT1", "OT2", "OC1", "OC2",
             "OXT"}


def _index(residue, atom) -> Optional[int]:
    for i, a in enumerate(residue.atoms):
        if a is atom:
            return i
    return None


def _heavy_neighbors(atom):
    return [nb for nb in getattr(atom, "bonds", []) if _element_of(nb) != "H"]


def _h_neighbors(atom):
    return [nb for nb in getattr(atom, "bonds", []) if _element_of(nb) == "H"]


def walk_reacting_atoms(residue, bond_atom_name: str,
                        boundary_depth: int = 2) -> Tuple[List[str], Optional[str]]:
    """Walk from the bond atom toward the backbone; return atoms to duplicate.

    Collects heavy sidechain atoms whose heavy-atom bond-distance from the bond
    atom is ``< boundary_depth`` (plus their attached hydrogens).  The first
    heavy atom at distance ``boundary_depth`` -- or the backbone -- is the
    *boundary*: it stays common and both duplicated branches reconnect to it.

    For a disulfide (bond atom ``SG``, depth 2) this yields reacting atoms
    ``SG, CB`` + their H (``HG1, HB1, HB2``) with boundary ``CA``.  For NOS-Lys
    (bond atom ``NZ``, depth 2): ``NZ, CE`` + H, boundary ``CD``.

    Parameters
    ----------
    residue : residue with ``.atoms`` carrying ``.bonds`` connectivity
    bond_atom_name : str
        Name of the atom that forms the inter-residue bond (e.g. ``'SG'``).
    boundary_depth : int
        Heavy-atom depth (from the bond atom) at which to stop duplicating.

    Returns
    -------
    (reacting_atom_names, boundary_atom_name)
    """
    by_name = {a.name: a for a in residue.atoms}
    if bond_atom_name not in by_name:
        raise KeyError(f"{bond_atom_name} not in residue {residue.resname}")
    start = by_name[bond_atom_name]

    # BFS over heavy atoms only, recording heavy-atom depth from the bond atom
    depth = {_index(residue, start): 0}
    order = [start]
    reacting_heavy = []
    boundary = None
    queue = [start]
    while queue:
        cur = queue.pop(0)
        d = depth[_index(residue, cur)]
        if cur.name in _BACKBONE:
            # reached the backbone: this is the boundary, do not duplicate it
            if boundary is None:
                boundary = cur.name
            continue
        if d >= boundary_depth:
            # this heavy atom is the boundary (stays common)
            if boundary is None:
                boundary = cur.name
            continue
        reacting_heavy.append(cur)
        for nb in _heavy_neighbors(cur):
            j = _index(residue, nb)
            if j not in depth:
                depth[j] = d + 1
                queue.append(nb)

    # reacting atoms = collected heavy atoms + their attached hydrogens
    reacting = []
    for hv in reacting_heavy:
        reacting.append(hv.name)
        for h in _h_neighbors(hv):
            reacting.append(h.name)

    return reacting, boundary


# ── split-atom (dual-topology) residue construction ───────────────────────────

def _split_name(name: str, branch: int) -> str:
    """Branch-tagged atom name: heavy 'CB' -> 'CB1'/'CB2'; H 'HB1' -> 'HB1'+tag.

    A per-branch suffix keeps names unique within the residue; the exact string
    is not physically meaningful (GROMACS only needs uniqueness + consistency
    between the structure and the RTP).
    """
    return f"{name}{branch}"


def _dual_atom(name, atomtype, q, m, x, symbol,
               nameB, atomtypeB, qB, mB) -> Atom:
    """Construct a pmx Atom carrying both alchemical states."""
    a = Atom(name=name, x=list(x), symbol=symbol)
    a.atomtype, a.q, a.m = atomtype, q, m
    a.nameB, a.atomtypeB, a.qB, a.mB = nameB, atomtypeB, qB, mB
    # mirror onto .type/.typeB too (apply_bond_morph's FF lookup reads these)
    a.type, a.typeB = atomtype, atomtypeB
    return a


def build_split_residue(wt_res, bound_res, bond_atom_name: str,
                        boundary_depth: int = 2) -> Tuple[List[Atom], List[Tuple[str, str]]]:
    """Auto-generate the split-atom dual-topology hybrid for a cross-link.

    Reusable for any cross-link (disulfide, NOS): duplicates the reacting
    sidechain (from :func:`walk_reacting_atoms`) into a free branch (WT params,
    real in A / dummy in B) and a bound branch (target params, dummy in A / real
    in B); non-reacting atoms are common morphs.

    Correspondence between *wt_res* and *bound_res* reacting atoms is by name.
    Atoms present only in the free residue (e.g. the thiol ``HG1``) live on the
    free branch only; atoms present only in the bound residue (e.g. a NOS
    bridging ``OD``) live on the bound branch only (appearing).

    Returns ``(atoms, bonds)`` where *atoms* are dual-state pmx Atoms and
    *bonds* are (name, name) pairs referencing them.
    """
    reacting, boundary = walk_reacting_atoms(wt_res, bond_atom_name, boundary_depth)
    reacting_set = set(reacting)
    wt_by = {a.name: a for a in wt_res.atoms}
    bound_by = {a.name: a for a in bound_res.atoms}

    atoms: List[Atom] = []

    # 1. common (non-reacting) atoms: morph WT -> bound (dummify B if absent)
    for a in wt_res.atoms:
        if a.name in reacting_set:
            continue
        b = bound_by.get(a.name)
        if b is not None:
            atoms.append(_dual_atom(a.name, a.atomtype, a.q, a.m, a.x, a.symbol,
                                    a.name, b.atomtype, b.q, b.m))
        else:  # WT-only common atom -> disappears
            atoms.append(_dual_atom(a.name, a.atomtype, a.q, a.m, a.x, a.symbol,
                                    a.name + ".gone", "DUM_" + a.atomtype,
                                    0.0, a.m))

    # 2. reacting atoms -> two decoupled branches
    for name in reacting:
        a = wt_by[name]
        n1 = _split_name(name, 1)
        # free branch: real in A, dummy in B
        atoms.append(_dual_atom(n1, a.atomtype, a.q, a.m, a.x, a.symbol,
                                n1 + ".gone", "DUM_" + a.atomtype, 0.0, a.m))
        b = bound_by.get(name)
        if b is not None:
            n2 = _split_name(name, 2)
            # bound branch: dummy in A, real in B.  It represents the SAME
            # physical atom as the free branch, so co-locate it on the WT atom's
            # coordinate (``a.x``) rather than the target's idealized geometry --
            # this is what the classic split placement (_set_conformation via
            # _SPLIT_ATOM_NAME_MAP) does, and it makes X1/X2 exactly coincident.
            atoms.append(_dual_atom(n2, "DUM_" + b.atomtype, 0.0, b.m, a.x,
                                    b.symbol, n2, b.atomtype, b.q, b.m))

    # 3. bound-only reacting atoms (e.g. NOS bridging OD): appearing on branch 2
    wt_names = {a.name for a in wt_res.atoms}
    for name, b in bound_by.items():
        if name in wt_names or name in reacting_set:
            continue
        n2 = _split_name(name, 2)
        atoms.append(_dual_atom(n2, "DUM_" + b.atomtype, 0.0, b.m, b.x, b.symbol,
                                n2, b.atomtype, b.q, b.m))

    bonds = _split_bonds(wt_res, bound_res, reacting_set, boundary)
    return atoms, bonds


class _SplitResidue:
    """Minimal residue view over split-atom dual-state atoms, for _make_mtp."""
    def __init__(self, atoms, resname, resnA, resnB):
        self.atoms = atoms
        self.resname = resname
        self.resnA = resnA
        self.resnB = resnB


class _RefAtom:
    """Placeholder for an inter-residue improper marker (``-C`` / ``+N``)."""
    def __init__(self, name):
        self.name = name


def _backbone_impropers(atom_map):
    """The standard protein backbone impropers, on the common backbone atoms.

    ``N -C CA HN`` and ``C CA +N O`` -- as in the shipped C2CD.  ``-C``/``+N``
    are inter-residue markers (placeholders); the rest resolve to this residue.
    """
    imps = []
    for quartet in (["N", "-C", "CA", "HN"], ["C", "CA", "+N", "O"]):
        own = [n for n in quartet if not n.startswith(("+", "-"))]
        if not all(n in atom_map for n in own):
            continue
        atoms = [atom_map[n] if n in atom_map else _RefAtom(n) for n in quartet]
        imps.append(atoms + ["default-A", "default-B"])
    return imps


def emit_crosslink_mtp(atoms: List[Atom], resname: str,
                       resnA: str, resnB: str) -> str:
    """Emit the split-atom cross-link residue as MTP text (reuses _make_mtp).

    Emits [ morphes ]/[ atoms ]/[ coords ] plus the common backbone impropers.
    The inter-residue bond is not part of the residue definition -- it is
    declared by the user at gentop (``pmx gentop --bond``).
    """
    from ..nsaa import _make_mtp
    for i, a in enumerate(atoms, 1):
        a.id = i
    r = _SplitResidue(atoms, resname, resnA, resnB)
    ii = _backbone_impropers({a.name: a for a in atoms})
    return "".join(_make_mtp(r, ii, [], rotations=None))


def emit_crosslink_rtp(atoms: List[Atom], bonds, resname: str,
                       resnA: str, resnB: str, ff: str) -> str:
    """Emit the split-atom cross-link residue as RTP text (reuses nsaa._make_rtp).

    *bonds* is the list of ``(name, name)`` intra-residue bond pairs from
    :func:`build_split_residue`; they are attached to the atoms as connectivity
    so ``_make_rtp`` can emit the ``[ bonds ]`` section.  The common backbone
    impropers, the peptide neighbour bond, and (for CHARMM) the backbone cmap
    are added so pdb2gmx accepts the residue.
    """
    from ..nsaa import _make_rtp
    by = {a.name: a for a in atoms}
    for a in atoms:
        a.bonds = []
    for n1, n2 in bonds:
        if n1 in by and n2 in by and by[n2] not in by[n1].bonds:
            by[n1].bonds.append(by[n2])
            by[n2].bonds.append(by[n1])
    for i, a in enumerate(atoms, 1):
        a.id = i
    r = _SplitResidue(atoms, resname, resnA, resnB)
    ii = _backbone_impropers(by)
    neigh = [("C", "+N")]                       # peptide bond to the next residue
    cmap = ["-C   N    CA     C    +N"] if "charmm" in ff else []
    return "".join(_make_rtp(r, ii, [], neigh, cmap, ff))


def decouple(wt_res, bound_res, bond_atom_name: str, boundary_depth: int = 2):
    """The unified 'decouple' operation: split a residue's reacting sidechain.

    Returns ``(atom_map, atoms, bonds)`` where *atom_map* is ``name -> Atom``
    over the split atoms (so callers reference the bound-branch ``*2`` atoms
    when building the cross-link bond-morph).  This is the single primitive the
    cross-link recipes below share.
    """
    atoms, bonds = build_split_residue(wt_res, bound_res, bond_atom_name,
                                       boundary_depth)
    return {a.name: a for a in atoms}, atoms, bonds


# ── model-level cross-link placement ──────────────────────────────────────────

def _superpose_by_name(target_res, ref_res) -> None:
    """Superpose *target_res* onto *ref_res* using shared heavy-atom names.

    Puts the bound target's coordinates in the WT residue's frame, so the split
    branches inherit sensible positions with no extra geometry step: the
    bound-branch reacting atoms land on their WT counterparts (a disulfide's
    ``CB2``/``SG2`` co-locate with the WT ``CB``/``SG``), and any bound-only atom
    (a NOS bridging ``OD``) is carried into the WT frame.  This is the same
    implicit-alignment the single-residue :func:`covalent.mutate` relies on --
    common atoms simply reuse the structure's existing positions.
    """
    from .mutate import _superpose_target
    ref_by = {a.name: i for i, a in enumerate(ref_res.atoms)}
    pairs = [(ref_by[a.name], j) for j, a in enumerate(target_res.atoms)
             if a.name in ref_by]
    _superpose_target(ref_res, target_res, pairs)


def _place_split_residue(model, wt_residue, split_atoms, hybrid_name):
    """Replace *wt_residue* in *model* with a residue built from *split_atoms*."""
    from ..molecule import Molecule
    new = Molecule()
    new.atoms = []
    new.resname = hybrid_name
    for a in split_atoms:
        a.resname = hybrid_name
        new.append(a)
    model.replace_residue(residue=wt_residue, new=new, bKeepResNum=True)
    return new




def decouple_mutate(model, resid, target, ff, react_atom, chain=None,
                    boundary_depth=2, inplace=False,
                    write_files=True, verbose=False):
    """Mutate one residue to a bound-state cross-link endpoint, with clawback.

    This is the per-residue half of a cross-link: it morphs *resid* to *target*
    (the reacted species, e.g. ``CN0``/``KN0``/``KS0``) AND splits the reacting
    sidechain (``react_atom``) into two branches so a forming inter-residue bond's
    1-4 exclusions are harmless.  It does *not* form the bond -- that is declared
    by the user later (``pmx gentop --bond``).  A SONOS bridge is just
    ``decouple_mutate`` on the Lys plus one on each Cys.

    Parameters
    ----------
    react_atom : str
        the reacting sidechain atom to split (``SG`` for Cys, ``NZ`` for Lys).

    Returns the generated hybrid residue name.

    Raises
    ------
    ValueError
        if the residue/atom selection is invalid, the residue is not a protein
        residue, or the cross-link target is not in the force field.
    """
    import warnings
    from copy import deepcopy
    from ..utils import list2file
    from ..library import pmx_data_file
    from ..nsaa import (_assign_rtp_entries, _get_residue_types, _ARN_CONTENT,
                        _write_atp_fnb, _assign_mass_atp)
    from .mutate import _load_rtp, _upsert_entry, _split_rtp_header
    from .target import TargetResidue

    rtp, ffpath = _load_rtp(ff)
    bCharmm = "charmm" in ff
    m = model if inplace else deepcopy(model)
    if m.unity == "nm":
        m.nm2a()

    # fetch_residue raises a clear ValueError for a bad/ambiguous resid or chain
    wt = m.fetch_residue(idx=resid, chain=chain)
    resnA = wt.resname

    # cross-links act on protein sidechains
    if getattr(wt, "moltype", None) not in (None, "", "protein", "unknown"):
        raise ValueError(
            "pmx mutate --react: residue %s %d is '%s' (%s), but a cross-link "
            "endpoint must be a protein residue." % (resnA, resid,
                                                     wt.moltype, resnA))

    _assign_rtp_entries(wt, rtp, None)
    if not any(getattr(a, "bonds", None) for a in wt.atoms):
        attach_bonds_from_rtp(wt, ff)

    # the reacting atom must exist in the residue being clawed back
    if not any(a.name == react_atom for a in wt.atoms):
        raise ValueError(
            "pmx mutate --react: residue %s %d has no atom '%s' to react at. "
            "Atoms present: %s" % (resnA, resid, react_atom,
                                   ", ".join(a.name for a in wt.atoms)))

    # a hydrogen-free WT leaves polar/backbone H as appearing dummies, which
    # pdb2gmx rejects (dangling improper) -- warn up front, as covalent.mutate does
    if not any(_element_of(a) == "H" for a in wt.atoms):
        warnings.warn(
            "pmx mutate --react: residue %s %d has no hydrogens -- the structure "
            "appears unprotonated; the hybrid's polar/backbone hydrogens become "
            "dummies and pdb2gmx will reject a dangling improper. Protonate the "
            "structure first (e.g. run gmx pdb2gmx on the WT and mutate that)."
            % (resnA, resid), stacklevel=2)

    # build the bound-state target; a missing residue is the endpoint-construction
    # gap (the bound species that defines the clawed-back branch).
    try:
        bound = TargetResidue.from_rtp(target, ff).residue
    except KeyError as e:
        raise ValueError(
            "pmx mutate --react: cross-link target '%s' is not in force field "
            "'%s'. The bound-state residue must exist before you can react into "
            "it -- e.g. CYS2/CYX for a disulfide. For NOS/SONOS/thioether, "
            "construct and ship the endpoint residue first (see CN0/KN0/KS0)."
            % (target, ff)) from e

    # the reacting atom must have a bound-state counterpart to split into
    if not any(a.name == react_atom for a in bound.atoms):
        raise ValueError(
            "pmx mutate --react: reacting atom '%s' is absent from the bound "
            "target '%s' (its atoms: %s) -- the clawback has no bound branch to "
            "split into." % (react_atom, target,
                             ", ".join(a.name for a in bound.atoms)))

    _assign_mass_atp(wt, bound, os.path.join(ffpath, "atomtypes.atp"))
    _superpose_by_name(bound, wt)
    amap, atoms, bonds = decouple(wt, bound, react_atom, boundary_depth)
    # Cross-link endpoint name = <one-letter WT code>2<target[:2]> (C2CD, C2CN,
    # K2KN, K2KS) -- matches residuetypes.dat, get_mtp_file routing and
    # library._nos_hybrids.  Use the one-letter CODE, not resnA[:1]: for LYS the
    # resname's first char is 'L' (-> wrong 'L2KN'), the code is 'K' (-> 'K2KN').
    letter = library._one_letter.get(resnA, resnA[:1])
    name = f"{letter}2{target[:2]}"[:4]

    if write_files:
        mtp_text = emit_crosslink_mtp(list(atoms), name, resnA, target)
        rtp_text = emit_crosslink_rtp(list(atoms), bonds, name, resnA, target, ff)
        _write_atp_fnb("atomtypes.atp", "ffnonbonded.itp",
                       _SplitResidue(atoms, name, resnA, target), ffpath)
        rtp_header, rtp_entry = _split_rtp_header(rtp_text, name)
        _upsert_entry(os.path.join(ffpath, "nsaa.rtp"), name, rtp_entry,
                      header=rtp_header)
        _upsert_entry(os.path.join(ffpath, "nsaa.mtp"), name, mtp_text)
        if bCharmm:
            arn = os.path.join(ffpath, "nsaa.arn")
            if not os.path.isfile(arn):
                list2file([_ARN_CONTENT], arn)
        fn_restypes = pmx_data_file("mutff/residuetypes.dat")
        if name not in _get_residue_types(fn_restypes):
            need_nl = (os.path.isfile(fn_restypes)
                       and os.path.getsize(fn_restypes) > 0
                       and open(fn_restypes, "rb").read()[-1:] != b"\n")
            with open(fn_restypes,
                      "a" if os.path.isfile(fn_restypes) else "w") as f:
                f.write(("\n" if need_nl else "") + f"{name:<6} Protein\n")

    _place_split_residue(m, wt, atoms, name)
    if verbose:
        print(f"log_> covalent.decouple_mutate: {resnA}{resid} -> {name} "
              f"(split at {react_atom}); MTP -> "
              f"{os.path.join(ffpath, 'nsaa.mtp')}")
    return name


def _split_bonds(wt_res, bound_res, reacting_set, boundary):
    """Reconnect the two branches: each reacting bond is duplicated per branch;
    a bond crossing the reacting/common boundary connects to both branches.

    Bound-only atoms (present only in the target, e.g. a NOS bridging ``OD``) are
    part of the bound branch too: they are branch-2 atoms, so their bonds are
    branch-2 bonds and their names carry the branch-2 suffix -- exactly like the
    reacting atoms.  (Without this, a bond such as the sulfenyl ``SG-OD`` would be
    emitted with the un-suffixed ``OD`` name, and an ``OD-HD`` bond between two
    bound-only atoms would be dropped entirely.)
    """
    out = []
    wt_names = {a.name for a in wt_res.atoms}
    # atoms living only on the bound branch: the reacting atoms (which split into
    # both branches) plus target-only atoms (which appear on branch 2 only).
    bound_only = {a.name for a in bound_res.atoms
                  if a.name not in wt_names and a.name not in reacting_set}
    branch2 = reacting_set | bound_only

    def side1(name):
        return _split_name(name, 1) if name in reacting_set else name

    def side2(name):
        return _split_name(name, 2) if name in branch2 else name

    seen = set()
    # WT bonds -> free branch (branch 1)
    for a in wt_res.atoms:
        for nb in getattr(a, "bonds", []):
            if a.name.startswith("+") or nb.name.startswith("+"):
                continue
            n1, n2 = sorted((a.name, nb.name))
            if (n1, n2, 1) in seen:
                continue
            seen.add((n1, n2, 1))
            out.append((side1(n1), side1(n2)))
    # bound bonds touching the bound branch (reacting or bound-only) -> branch 2
    seen2 = set()
    for a in bound_res.atoms:
        for nb in getattr(a, "bonds", []):
            if a.name.startswith("+") or nb.name.startswith("+"):
                continue
            if a.name not in branch2 and nb.name not in branch2:
                continue
            n1, n2 = sorted((a.name, nb.name))
            if (n1, n2) in seen2:
                continue
            seen2.add((n1, n2))
            out.append((side2(n1), side2(n2)))
    return out
