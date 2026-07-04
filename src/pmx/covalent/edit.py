"""Inter-residue covalent edits: a declarative unification of the bespoke
disulfide / NOS / terminal-deletion bonded-term logic in :mod:`pmx.alchemy`.

A covalent modification that spans more than one residue (an alchemical
disulfide, an N-O-S bridge, a terminal deletion) is represented as a
:class:`CovalentEdit`: a set of participating residues plus a list of
:class:`InterResidueBondMorph` -- cross-residue bonded interactions
(bonds / angles / dihedrals) that are toggled between the A and B states.

Each morph carries a ``mode``:

* ``form``  -- A-state off (k=0), B-state on   (e.g. disulfide/NOS formation)
* ``break`` -- A-state on,       B-state off   (bond breaking)
* ``perma`` -- both states on
* ``decouple`` -- zero the B-state of *existing* dihedrals spanning a boundary
  (terminal deletion: the target residue is fully dummy in state B)

Force constants come from the force field at apply time (via
``topol.BondedParams``), exactly as the original routines did -- the morph
only declares *which* atoms and interactions participate, not the numbers.
This is what the ``[ inter_bonds ]`` MTP section will serialize.
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


@dataclass
class CovalentEdit:
    """A covalent modification spanning one or more residues.

    The unified surface: compose a modification from the operations
    :meth:`decouple` (split a reacting sidechain for a forming bond) and
    :meth:`bond` (toggle an inter-residue bond).  A single-residue mutation is a
    :class:`CovalentEdit` with neither.

    ``residues`` are the participating residue names; ``bond_morphs`` the
    inter-residue bonded transformations to apply; ``decoupled`` the split-atom
    residues produced by :meth:`decouple` (each ``(resname, atoms, bonds)``).
    """
    residues: List[str] = field(default_factory=list)
    bond_morphs: List["InterResidueBondMorph"] = field(default_factory=list)
    decoupled: List = field(default_factory=list)

    def decouple(self, wt_res, bound_res, bond_atom_name, boundary_depth=2):
        """Split *wt_res*'s reacting sidechain (WT->bound). Returns name->Atom.

        The returned map lets a caller reference the bound-branch ``*2`` atoms
        when building the cross-link bond-morph.
        """
        from .crosslink import decouple as _decouple
        amap, atoms, bonds = _decouple(wt_res, bound_res, bond_atom_name,
                                       boundary_depth)
        self.residues.append(wt_res.resname)
        self.decoupled.append((wt_res.resname, atoms, bonds))
        return amap

    def bond(self, morph: "InterResidueBondMorph") -> "CovalentEdit":
        """Register an inter-residue bond-morph; returns self for chaining."""
        self.bond_morphs.append(morph)
        return self


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
    for it in morph.interactions:
        types = [_state_type(a, morph.mode) for a in it.atoms]
        if it.kind == "bond":
            func, x0, k = bp.get_bond_param(*types)
            topol.bonds.append(list(it.atoms) +
                               [func, [func, x0, k * fA], [func, x0, k * fB]])
        elif it.kind == "angle":
            parm = bp.get_angle_param(*types)
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
            for func, phi0, kphi, mult in bp.get_dihedral_param(*types, 9):
                topol.dihedrals.append(list(it.atoms) + [
                    func,
                    [func, phi0, kphi * fA, mult],
                    [func, phi0, kphi * fB, mult],
                ])
                _zero_B_if_has_dummy_dihedral(topol.dihedrals[-1])
        else:
            raise ValueError(f"unknown interaction kind: {it.kind!r}")


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


def apply_edit(topol, edit: CovalentEdit) -> None:
    """Apply all bond morphs of a covalent edit to *topol*."""
    for morph in edit.bond_morphs:
        apply_bond_morph(topol, morph)


# ── constructors (declarative equivalents of the alchemy routines) ────────────

def disulfide_morph(cys1: Dict, cys2: Dict, mode: str = "form",
                    ca="CA", cb="CB", sg="SG") -> InterResidueBondMorph:
    """SG-SG bond + CB-SG-SG angles + SS dihedrals (cf. _add_ss_bonded_terms).

    Each of *cys1*, *cys2* is a name->Atom dict.  The atom-name keys default to
    the graph-engine convention (``CA``/``CB``/``SG`` -- a plain CYS->CYS2 morph
    through ``build_hybrid_residue``); pass ``cb='CB2'``/``sg='SG2'`` for the
    legacy split-atom encoding.
    """
    I = Interaction
    ca1, cb1, sg1 = cys1[ca], cys1[cb], cys1[sg]
    ca2, cb2, sg2 = cys2[ca], cys2[cb], cys2[sg]
    return InterResidueBondMorph(mode, [
        I("bond", [sg1, sg2]),
        I("angle", [cb1, sg1, sg2]),
        I("angle", [cb2, sg2, sg1]),
        I("dihedral", [cb1, sg1, sg2, cb2]),
        I("dihedral", [ca1, cb1, sg1, sg2]),
        I("dihedral", [ca2, cb2, sg2, sg1]),
    ])


def nos_morph(cys: Dict, lys: Dict, mode: str = "form") -> InterResidueBondMorph:
    """NZ-OD2 bond + angles + dihedrals of an N-O-S bridge (cf. _add_nos_terms).

    *cys* keys: CA, CB2, SG2, OD2.  *lys* keys: CD, CE2, NZ2, HZ21, HE21, HE22.
    Optional atoms (HZ21/HE21/HE22) are included when present.
    """
    I = Interaction
    nz, od2, sg2 = lys["NZ2"], cys["OD2"], cys["SG2"]
    ce, cd, cb2 = lys["CE2"], lys["CD"], cys["CB2"]
    hz1, he1, he2 = lys.get("HZ21"), lys.get("HE21"), lys.get("HE22")

    interactions = [I("bond", [nz, od2])]
    interactions.append(I("angle", [ce, nz, od2]))
    if hz1 is not None:
        interactions.append(I("angle", [hz1, nz, od2]))
    interactions.append(I("angle", [nz, od2, sg2]))
    interactions += [
        I("dihedral", [cd, ce, nz, od2]),
        I("dihedral", [ce, nz, od2, sg2]),
        I("dihedral", [nz, od2, sg2, cb2]),
    ]
    if hz1 is not None:
        interactions.append(I("dihedral", [hz1, nz, od2, sg2]))
    if he1 is not None:
        interactions.append(I("dihedral", [he1, ce, nz, od2]))
    if he2 is not None:
        interactions.append(I("dihedral", [he2, ce, nz, od2]))
    return InterResidueBondMorph(mode, interactions)
