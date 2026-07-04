"""Topology-first atom mapping for covalent residue mutations.

The classic pmx ligand mapping (``pmx.ligand_alchemy.LigandAtomMapping`` /
``pmx.nsaa.domapping``) is built for *arbitrary* ligand pairs with no atom-name
correspondence, so it relies on a 3-D maximum-common-substructure search with
distance-based pruning.  For a residue -> RTP-target mutation, however, both
residues come from the *same* force field's naming convention, so atom names
(``N``, ``CA``, ``CB``, ``OG``, ``HN`` ...) are an exact, free correspondence.

:func:`map_residues` uses that directly, in tiers:

1. **Name + element anchor** -- atoms sharing a (canonicalized name, element)
   are the same atom.  For same-force-field residues this deterministically
   covers the whole shared core, hydrogens included, with no 3-D dependence.
2. **Graph-MCS fallback** (planned) -- for atoms whose names do *not*
   correspond (novel groups, cross-force-field, arbitrary ligand attachments),
   a 2-D graph MCS over the residual atoms, seeded by the Tier-1 anchor.

Unmatched atoms fall out as dummies: WT-only -> *disappearing* (real in state
A, dummy in B), target-only -> *appearing* (dummy in A, real in B).  Because
correspondence is by identity rather than a post-hoc patch, polar hydrogens
(e.g. ``HN``) map naturally and there is no silent length-mismatch failure.

Coordinates are intentionally *not* used for the mapping; they are only needed
downstream to place the dummy atoms.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


# ── name / element helpers ────────────────────────────────────────────────────

def _canon_name(name: str) -> str:
    """Canonicalize an atom name across PDB/GROMACS hydrogen conventions.

    A leading digit (PDB style, e.g. ``1HB``) is rotated to the end to match the
    GROMACS style (``HB1``).  Non-hydrogen names are unaffected in practice.
    """
    n = name.strip()
    if len(n) > 1 and n[0].isdigit():
        return n[1:] + n[0]
    return n


# element symbol -> canonical form; only the leading letter is used for the
# name-based fallback, which is unambiguous for standard protein/nucleic atoms
# (Cl/Br etc. do not occur as residue atom names).
def _element_of(atom) -> str:
    """Best-effort chemical element for a residue atom."""
    sym = (getattr(atom, "symbol", "") or "").strip()
    if sym:
        return sym[0].upper() + sym[1:].lower()
    n = atom.name.strip()
    if n and n[0].isdigit():
        n = n[1:]
    return n[0].upper() if n else "C"


# ── result container ──────────────────────────────────────────────────────────

@dataclass
class ResidueMapping:
    """Correspondence between a WT residue (state A) and a target (state B).

    ``pairs`` are ``(i, j)`` index pairs into ``res1.atoms`` / ``res2.atoms``.
    ``disappear`` are res1 indices with no partner (real in A, dummy in B);
    ``appear`` are res2 indices with no partner (dummy in A, real in B).
    ``method`` records how each res1 atom was matched (``'name'`` / ``'mcs'``).
    """

    res1: object
    res2: object
    pairs: List[Tuple[int, int]] = field(default_factory=list)
    disappear: List[int] = field(default_factory=list)
    appear: List[int] = field(default_factory=list)
    method: Dict[int, str] = field(default_factory=dict)

    def pairs1(self) -> List[int]:
        return [i for i, _ in self.pairs]

    def pairs2(self) -> List[int]:
        return [j for _, j in self.pairs]

    def named_pairs(self) -> List[Tuple[str, str]]:
        return [(self.res1.atoms[i].name, self.res2.atoms[j].name)
                for i, j in self.pairs]

    def summary(self) -> str:
        return (f"{len(self.pairs)} mapped, "
                f"{len(self.disappear)} disappear, {len(self.appear)} appear")


# ── mapping ───────────────────────────────────────────────────────────────────

def map_residues(res1, res2, element_guard: bool = True,
                 mcs_fallback: bool = False) -> ResidueMapping:
    """Map two residues by atom name+element, classifying the remainder as dummies.

    Parameters
    ----------
    res1, res2 : residue-like
        Objects with an ``.atoms`` list of pmx :class:`~pmx.atom.Atom`.
        ``res1`` is the WT (state A), ``res2`` the target (state B).
    element_guard : bool
        Require matched names to share an element (guards against rare
        cross-element name collisions).  Default True.
    mcs_fallback : bool
        Run a graph-MCS pass over atoms left unmatched by the name anchor.
        Not implemented yet; see module docstring.

    Returns
    -------
    ResidueMapping
    """
    # index res2 atoms by canonical name (preserve order for determinism)
    r2_by_canon: Dict[str, List[int]] = {}
    for j, a in enumerate(res2.atoms):
        r2_by_canon.setdefault(_canon_name(a.name), []).append(j)

    used2: set = set()
    pairs: List[Tuple[int, int]] = []
    method: Dict[int, str] = {}
    unmatched1: List[int] = []

    # Tier 1: name + element anchor
    for i, a in enumerate(res1.atoms):
        canon = _canon_name(a.name)
        chosen = None
        for j in r2_by_canon.get(canon, []):
            if j in used2:
                continue
            if element_guard and _element_of(a) != _element_of(res2.atoms[j]):
                continue
            chosen = j
            break
        if chosen is None:
            unmatched1.append(i)
        else:
            pairs.append((i, chosen))
            used2.add(chosen)
            method[i] = "name"

    unmatched2 = [j for j in range(len(res2.atoms)) if j not in used2]

    # Tier 2: graph-MCS fallback over the residual atoms (planned)
    if mcs_fallback and unmatched1 and unmatched2:
        raise NotImplementedError(
            "graph-MCS fallback (Tier 2) is not implemented yet; "
            "call with mcs_fallback=False for the name-anchor mapping.")

    return ResidueMapping(res1=res1, res2=res2, pairs=pairs,
                          disappear=unmatched1, appear=unmatched2,
                          method=method)
