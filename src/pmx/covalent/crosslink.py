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

from typing import List, Optional, Tuple

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
