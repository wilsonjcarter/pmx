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

from typing import Dict, List, Optional, Tuple

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
            # bound branch: dummy in A, real in B
            atoms.append(_dual_atom(n2, "DUM_" + b.atomtype, 0.0, b.m, b.x,
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
                       resnA: str, resnB: str, inter_bonds=None) -> str:
    """Emit the split-atom cross-link residue as MTP text (reuses _make_mtp).

    Emits [ morphes ]/[ atoms ]/[ coords ] plus the common backbone impropers,
    and -- when *inter_bonds* is given -- an ``[ inter_bonds ]`` section listing
    the cross-residue bonds this residue forms as ``(mode, self, partner)``
    tuples (``gen_hybrid_top`` derives the induced angles/dihedrals at apply
    time).
    """
    from ..nsaa import _make_mtp
    for i, a in enumerate(atoms, 1):
        a.id = i
    r = _SplitResidue(atoms, resname, resnA, resnB)
    ii = _backbone_impropers({a.name: a for a in atoms})
    return "".join(_make_mtp(r, ii, [], rotations=None, inter_bonds=inter_bonds))


def build_crosslink_mtp(wt_res, bound_res, bond_atom_name: str,
                        resname: str, boundary_depth: int = 2,
                        inter_bonds=None) -> str:
    """Walk + duplicate + emit: the split-atom cross-link hybrid as MTP text."""
    atoms, _bonds = build_split_residue(wt_res, bound_res, bond_atom_name,
                                        boundary_depth)
    return emit_crosslink_mtp(atoms, resname, wt_res.resname, bound_res.resname,
                              inter_bonds=inter_bonds)


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


def disulfide_edit(cys1_res, cys2_res, bound1_res, bound2_res, mode="form"):
    """Alchemical disulfide = decouple x2 + SG-SG bond on the bound branch.

    ``cysN_res`` = WT cysteines; ``boundN_res`` = disulfide-Cys targets
    (``CYS2``/``CYX``).  Returns ``(edit, split1, split2)``.
    """
    from .edit import CovalentEdit, disulfide_morph
    edit = CovalentEdit()
    a1 = edit.decouple(cys1_res, bound1_res, "SG")
    a2 = edit.decouple(cys2_res, bound2_res, "SG")
    edit.bond(disulfide_morph(a1, a2, mode=mode, cb="CB2", sg="SG2"))
    return edit, edit.decoupled[0][1:], edit.decoupled[1][1:]


def nos_edit(cys_res, lys_res, nos_cys_res, nos_lys_res, mode="form"):
    """Alchemical N-O-S bridge = decouple Cys + decouple Lys + NZ-OD/OD-SG bonds.

    Same two primitives as the disulfide (``decouple`` + a bond-morph), just a
    different recipe: the Cys reacts at ``SG`` (its bound target gains the
    bridging ``OD``), the Lys reacts at ``NZ``.  ``nos_morph`` already reads the
    bound-branch ``*2`` names, so no remapping is needed.  ``nos_cys_res`` /
    ``nos_lys_res`` are the bound NOS targets (e.g. ``CSO``-derived for Cys).
    """
    from .edit import CovalentEdit, nos_morph
    edit = CovalentEdit()
    ac = edit.decouple(cys_res, nos_cys_res, "SG")
    al = edit.decouple(lys_res, nos_lys_res, "NZ")
    edit.bond(nos_morph(ac, al, mode=mode))
    return edit, edit.decoupled[0][1:], edit.decoupled[1][1:]


def _split_bonds(wt_res, bound_res, reacting_set, boundary):
    """Reconnect the two branches: each reacting bond is duplicated per branch;
    a bond crossing the reacting/common boundary connects to both branches."""
    out = []

    def side(name, branch):
        return _split_name(name, branch) if name in reacting_set else name

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
            out.append((side(n1, 1), side(n2, 1)))
    # bound bonds among reacting/boundary atoms -> bound branch (branch 2)
    seen2 = set()
    for a in bound_res.atoms:
        for nb in getattr(a, "bonds", []):
            if a.name.startswith("+") or nb.name.startswith("+"):
                continue
            # only add branch-2 duplicates for bonds touching the reacting set
            if a.name not in reacting_set and nb.name not in reacting_set:
                continue
            n1, n2 = sorted((a.name, nb.name))
            if (n1, n2) in seen2:
                continue
            seen2.add((n1, n2))
            out.append((side(n1, 2), side(n2, 2)))
    return out
