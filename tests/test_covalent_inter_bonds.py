#!/usr/bin/env python
"""Tests for inter-residue cross-link bonds at gentop.

Inter-residue bonds are declared by the user at gentop (``pmx gentop --bond
r1:atom-r2:atom``) -- there is no ``[ inter_bonds ]`` MTP section and no resname
auto-detection.  ``_apply_inter_bonds`` resolves each pair against the parsed
topology (plain atom name -> bound branch), then derives the induced
angles/dihedrals from the union bond graph.  These tests use hand-built objects
(no RDKit / force field) so they pin the apply logic directly.
"""
from pmx.atom import Atom


def _atom(name, i, t="CT2", tB="DUM_CT2"):
    a = Atom(name=name, id=i)
    a.atomtype, a.atomtypeB = t, tB
    a.nameB, a.q, a.qB, a.m, a.mB = name, 0.0, 0.0, 12.0, 12.0
    a.x = [0.0, 0.0, 0.0]
    a.type, a.typeB = t, tB
    return a


class _R:
    def __init__(self, atoms, resname, resnA, resnB):
        self.atoms, self.resname, self.resnA, self.resnB = atoms, resname, resnA, resnB


def test_writer_omits_inter_bonds_section():
    """The MTP writer never emits an [ inter_bonds ] section -- inter-residue
    bonds live at gentop (--bond), not in the residue definition."""
    from pmx.nsaa import _make_mtp
    r = _R([_atom("SG2", 1)], "C2CY", "CYS", "CYS2")
    txt = "".join(_make_mtp(r, [], []))
    assert "[ inter_bonds ]" not in txt


# ── consumer: _apply_inter_bonds over a synthetic topology ────────────────────

class _BP:
    def get_bond_param(self, *t):     return (1, 2.05, 300.0)
    def get_angle_param(self, *t):    return (5, 110.0, 50.0, 0.0, 0.0)
    def get_dihedral_param(self, *t): return [(9, 180.0, 2.0, 2)]


class _Topol:
    def __init__(self, bonds):
        self.bonds = [list(b) for b in bonds]
        self.angles, self.dihedrals = [], []
        self.BondedParams = _BP()


class _Res:
    def __init__(self, resname, atoms, resid=0):
        self.resname, self.atoms, self.id = resname, atoms, resid


def _mkres(resname, base_id, resid):
    names = ["CA", "CB2", "SG2", "HB12", "HB22"]
    atoms = {}
    for k, n in enumerate(names):
        a = Atom(name=n, id=base_id + k)
        a.type = a.typeB = "X"
        a.resname = resname
        atoms[n] = a
    return _Res(resname, list(atoms.values()), resid), atoms


def test_gentop_bond_forms_disulfide():
    """gentop --bond SG-SG forms one dual-state bond + graph-induced terms.

    The plain atom name 'SG' maps to the bound branch 'SG2'."""
    from pmx.alchemy import _apply_inter_bonds
    r1, a1 = _mkres("C2QN", 1, 32)
    r2, a2 = _mkres("C2QN", 100, 35)
    bonds = []
    for a in (a1, a2):
        bonds += [[a["CA"], a["CB2"], 1], [a["CB2"], a["SG2"], 1],
                  [a["CB2"], a["HB12"], 1], [a["CB2"], a["HB22"], 1]]
    topol = _Topol(bonds)
    _apply_inter_bonds(topol, [r1, r2], [
        {"resid_a": 32, "atom_a": "SG", "resid_b": 35, "atom_b": "SG",
         "mode": "form"}])
    assert len(topol.bonds) == len(bonds) + 1          # one SG2-SG2' bond
    nb = topol.bonds[-1]
    assert {nb[0].name, nb[1].name} == {"SG2"}
    assert nb[3][2] == 0.0 and nb[4][2] == 300.0       # state A off, B on
    assert len(topol.angles) == 2
    assert len(topol.dihedrals) == 7                   # graph-complete


def test_gentop_bond_multibond_sonos():
    """Two --bond meeting at the shared Lys NZ2 (SONOS): both bonds plus the
    shared-nitrogen ODa-NZ-ODb cross-angle; partners disambiguated by resid."""
    from pmx.alchemy import _apply_inter_bonds

    def lys(resid):
        d = {n: Atom(name=n, id=resid * 10 + k)
             for k, n in enumerate(["CE2", "NZ2", "HZ21"])}
        for a in d.values():
            a.type = a.typeB = "X"; a.resname = "L2ZS"
        return _Res("L2ZS", list(d.values()), resid), d

    def cys(resid):
        d = {n: Atom(name=n, id=resid * 10 + k)
             for k, n in enumerate(["CB2", "SG2", "OD2"])}
        for a in d.values():
            a.type = a.typeB = "X"; a.resname = "C2QN"
        return _Res("C2QN", list(d.values()), resid), d

    rl, l = lys(72); rca, ca = cys(73); rcb, cb = cys(90)
    bonds = [[l["CE2"], l["NZ2"], 1], [l["NZ2"], l["HZ21"], 1],
             [ca["CB2"], ca["SG2"], 1], [ca["SG2"], ca["OD2"], 1],
             [cb["CB2"], cb["SG2"], 1], [cb["SG2"], cb["OD2"], 1]]
    topol = _Topol(bonds)
    _apply_inter_bonds(topol, [rl, rca, rcb], [
        {"resid_a": 72, "atom_a": "NZ", "resid_b": 73, "atom_b": "OD",
         "mode": "form"},
        {"resid_a": 72, "atom_a": "NZ", "resid_b": 90, "atom_b": "OD",
         "mode": "form"}])
    assert len(topol.bonds) == len(bonds) + 2          # two forming bonds
    cross = [ang for ang in topol.angles
             if ang[1].name == "NZ2"
             and ang[0].name == "OD2" and ang[2].name == "OD2"]
    assert cross                                       # ODa-NZ2-ODb cross-angle


def test_gentop_bond_noop_without_bonds():
    from pmx.alchemy import _apply_inter_bonds
    r1, a1 = _mkres("SOME", 1, 5)
    topol = _Topol([[a1["CA"], a1["CB2"], 1]])
    _apply_inter_bonds(topol, [r1], None)              # no --bond -> no-op
    assert topol.angles == [] and topol.dihedrals == [] and len(topol.bonds) == 1


def test_apply_bond_morph_skips_missing_param_with_warning():
    """A generated induced term with no FF parameter is skipped + warned, not a
    crash (surfaces gaps cleanly when params are incomplete/placeholder)."""
    import warnings
    from pmx.covalent.edit import (apply_bond_morph, InterResidueBondMorph,
                                   Interaction)

    class _BPNone:
        def get_bond_param(self, *t):     return (1, 2.0, 300.0)
        def get_angle_param(self, *t):    return None      # missing
        def get_dihedral_param(self, *t): return []

    class _T:
        def __init__(self):
            self.bonds, self.angles, self.dihedrals = [], [], []
            self.BondedParams = _BPNone()

    a = _atom("NZ2", 1, "NG301", "NG301")
    b = _atom("OD2", 2, "OG301", "OG301")
    c = _atom("OD2", 3, "OG301", "OG301")
    morph = InterResidueBondMorph("form", [
        Interaction("bond", [a, b]),
        Interaction("angle", [b, a, c]),          # OD-NZ-OD, no FF param
    ])
    topol = _T()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        apply_bond_morph(topol, morph)
    assert len(topol.bonds) == 1                   # the bond still formed
    assert len(topol.angles) == 0                  # the unparameterised angle skipped
    assert any("SKIPPED" in str(x.message) for x in w)
