"""Inter-residue covalent edits: the declarative bond-morph applied when a
cross-link forms.

A cross-residue bond (an alchemical disulfide, an N-O-S bridge, a terminal
deletion) is represented as an :class:`InterResidueBondMorph` -- the bonded
interactions (bond / angles / dihedrals) that are toggled between the A and B
states.  The forming bond is declared by the user (``pmx gentop --bond``) and
the induced angles/dihedrals are derived from the molecular graph by
:func:`induced_bond_morph`.

Each morph carries a ``mode``:

* ``form``  -- A-state off (k=0), B-state on   (e.g. disulfide/NOS formation)
* ``break`` -- A-state on,       B-state off   (bond breaking)
* ``perma`` -- both states on
* ``decouple`` -- zero the B-state of *existing* dihedrals spanning a boundary
  (terminal deletion: the target residue is fully dummy in state B)

Force constants come from the force field at apply time (via
``topol.BondedParams``), exactly as the original routines did -- the morph
only declares *which* atoms and interactions participate, not the numbers.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List

# Reuse the dummy-zeroing helpers from alchemy (behaviour must stay identical).
from ..alchemy import _zero_B_if_has_dummy_angle, _zero_B_if_has_dummy_dihedral


# mode -> (A-state force-constant factor, B-state factor)
_MODE_FACTORS = {
    "form":  (0.0, 1.0),
    "break": (1.0, 0.0),
    "perma": (1.0, 1.0),
}


@dataclass
class Interaction:
    """One cross-residue bonded interaction over resolved Atom objects."""
    kind: str            # 'bond' | 'angle' | 'dihedral'
    atoms: List          # list of pmx Atom objects


@dataclass
class InterResidueBondMorph:
    """A set of inter-residue interactions toggled between states by *mode*."""
    mode: str                     # 'form' | 'break' | 'perma' | 'decouple'
    interactions: List[Interaction] = field(default_factory=list)


# ── applying a morph to a topology ────────────────────────────────────────────

def _state_type(atom, mode):
    """Which atom-type drives the FF lookup: B-state for 'form', else A-state."""
    return atom.typeB if mode == "form" else atom.type


def apply_bond_morph(topol, morph: InterResidueBondMorph) -> None:
    """Append/modify *topol* bonded terms for one inter-residue morph.

    Reproduces the behaviour of alchemy._add_ss_bonded_terms /
    _add_nos_terms (for the additive modes) declaratively.
    """
    if morph.mode == "decouple":
        _apply_decouple(topol, morph)
        return

    fA, fB = _MODE_FACTORS[morph.mode]
    bp = topol.BondedParams
    missing = []
    for it in morph.interactions:
        types = [_state_type(a, morph.mode) for a in it.atoms]
        if it.kind == "bond":
            parm = bp.get_bond_param(*types)
            if parm is None:
                missing.append(("bond", tuple(types))); continue
            func, x0, k = parm
            topol.bonds.append(list(it.atoms) +
                               [func, [func, x0, k * fA], [func, x0, k * fB]])
        elif it.kind == "angle":
            parm = bp.get_angle_param(*types)
            if parm is None:
                missing.append(("angle", tuple(types))); continue
            if len(parm) != 5:
                parm = list(parm) + [0, 0]
            func, th0, kth, r0_ub, kub = parm
            topol.angles.append(list(it.atoms) + [
                func,
                [func, th0, kth * fA, r0_ub, kub * fA],
                [func, th0, kth * fB, r0_ub, kub * fB],
            ])
            _zero_B_if_has_dummy_angle(topol.angles[-1])
        elif it.kind == "dihedral":
            dparm = bp.get_dihedral_param(*types, 9)
            if not dparm:
                missing.append(("dihedral", tuple(types))); continue
            for func, phi0, kphi, mult in dparm:
                topol.dihedrals.append(list(it.atoms) + [
                    func,
                    [func, phi0, kphi * fA, mult],
                    [func, phi0, kphi * fB, mult],
                ])
                _zero_B_if_has_dummy_dihedral(topol.dihedrals[-1])
        else:
            raise ValueError(f"unknown interaction kind: {it.kind!r}")
    if missing:
        import warnings
        uniq = sorted(set(missing))
        warnings.warn(
            "apply_bond_morph: %d induced term(s) have no force-field parameter "
            "and were SKIPPED (the cross-link topology is incomplete until these "
            "are added): %s" % (len(uniq),
                                 "; ".join("%s %s" % (k, "-".join(t))
                                           for k, t in uniq)), stacklevel=2)


def _apply_decouple(topol, morph: InterResidueBondMorph) -> None:
    """Zero the B-state of existing dihedrals listed in the morph (terminal)."""
    target_ids = {frozenset(a.id for a in it.atoms)
                  for it in morph.interactions if it.kind == "dihedral"}
    for dih in topol.dihedrals:
        if frozenset({dih[0].id, dih[1].id, dih[2].id, dih[3].id}) not in target_ids:
            continue
        if dih[5] == "NULL" or dih[6] == "NULL":
            continue
        func = dih[4]
        if func == 3:
            dih[6] = [func, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        elif func == 2:
            dih[6] = [func, dih[6][1], 0.0]
        else:
            dih[6] = [func, dih[6][1], 0.0, dih[6][-1]]


# ── the generic, graph-derived bond-morph (the unified core) ──────────────────

def induced_bond_morph(new_bonds, adjacency: Dict, mode: str = "form"
                       ) -> InterResidueBondMorph:
    """Build an inter-residue bond-morph from *only* the forming bond(s).

    A cross-link is fully specified by which bonds form: the induced angles and
    dihedrals follow deterministically from the molecular graph -- exactly as
    GROMACS generates bonded terms from the bond list, only here in dual (A/B)
    state form.

    Given the forming bond(s) and the *union* adjacency of the participating
    residues (which must already include *new_bonds*), emit:

    * each new bond;
    * every angle ``x-a-b`` / ``a-b-y`` incident on a new edge ``a-b``;
    * every proper dihedral whose 4-atom path traverses a new edge.

    Terms are canonicalised so each physical interaction appears once, and their
    order is deterministic (graph order of *adjacency*), so the output is
    reproducible.  Force constants are resolved later, in :func:`apply_bond_morph`.

    Parameters
    ----------
    new_bonds : iterable of (Atom, Atom)
        The inter-residue bonds that form (e.g. ``[(SG, SG')]`` for a disulfide,
        ``[(NZ, OD)]`` for a NOS bridge).
    adjacency : dict {Atom: set(Atom)}
        The union bond graph of the participating residues, *including*
        *new_bonds*.
    mode : str
        ``'form'`` | ``'break'`` | ``'perma'``.
    """
    new_edges = {frozenset((a, b)) for a, b in new_bonds}

    def _order(node):
        # deterministic neighbour ordering (by atom id, falling back to name)
        return sorted(adjacency.get(node, ()),
                      key=lambda n: (getattr(n, "id", 0) or 0, n.name))

    def _traverses_new(path):
        return any(frozenset(pair) in new_edges
                   for pair in zip(path, path[1:]))

    interactions, seen_ang, seen_dih = [], set(), set()

    # bonds
    for a, b in new_bonds:
        interactions.append(Interaction("bond", [a, b]))

    # angles + dihedrals: enumerate simple paths of length 2 (angles) and 3
    # (dihedrals) hanging off each new edge, keeping only those that traverse it.
    for a, b in new_bonds:
        for mid_lo, mid_hi in ((a, b), (b, a)):
            # angle  x - mid_lo - mid_hi
            for x in _order(mid_lo):
                if x is mid_hi:
                    continue
                key = (frozenset((x, mid_hi)), mid_lo)
                if key in seen_ang:
                    continue
                seen_ang.add(key)
                interactions.append(Interaction("angle", [x, mid_lo, mid_hi]))

        # dihedrals: every 4-atom simple path containing edge (a, b)
        for x in _order(a):                       # x - a - b - y  (edge central)
            if x is b:
                continue
            for y in _order(b):
                if y is a or y is x:
                    continue
                _add_dih(interactions, seen_dih, [x, a, b, y])
        for x in _order(a):                       # y - x - a - b  (edge terminal)
            if x is b:
                continue
            for y in _order(x):
                if y is a or y is b:
                    continue
                _add_dih(interactions, seen_dih, [y, x, a, b])
        for x in _order(b):                       # a - b - x - y  (edge terminal)
            if x is a:
                continue
            for y in _order(x):
                if y is b or y is a:
                    continue
                _add_dih(interactions, seen_dih, [a, b, x, y])

    return InterResidueBondMorph(mode, interactions)


def _add_dih(interactions, seen, quartet):
    """Append a proper dihedral once (canonicalised against its reverse)."""
    ids = tuple(a.id for a in quartet)
    key = min(ids, tuple(reversed(ids)))
    if key in seen:
        return
    seen.add(key)
    interactions.append(Interaction("dihedral", list(quartet)))
