#!/usr/bin/env python
"""Tests for pmx.covalent.edit -- the declarative inter-residue bond-morph IR
(no RDKit or force field needed -- uses a mock topology).
"""
import pytest

from pmx.covalent.edit import (disulfide_morph, nos_morph, apply_bond_morph,
                               InterResidueBondMorph, CovalentEdit)


class _BP:
    def get_bond_param(self, *t):     return (1, 2.05, 100.0)
    def get_angle_param(self, *t):    return (5, 110.0, 50.0, 0.0, 0.0)
    def get_dihedral_param(self, *t): return [(9, 180.0, 2.0, 2), (9, 0.0, 1.0, 3)]


class _Topol:
    def __init__(self):
        self.bonds, self.angles, self.dihedrals = [], [], []
        self.BondedParams = _BP()


class _At:
    def __init__(self, name, i):
        self.name, self.id = name, i
        self.type = self.typeB = "SS"   # no 'DUM' -> no B-state zeroing


def _cys(base):
    return {"CA": _At("CA", base), "CB2": _At("CB2", base + 1),
            "SG2": _At("SG2", base + 2)}


def _sig(entry, n_atom):
    return (tuple(a.name for a in entry[:n_atom]), entry[n_atom],
            tuple(entry[n_atom + 1]), tuple(entry[n_atom + 2]))


def _dump(topol):
    return ([_sig(b, 2) for b in topol.bonds]
            + [_sig(a, 3) for a in topol.angles]
            + [_sig(d, 4) for d in topol.dihedrals])


@pytest.mark.parametrize("mode", ["form", "break", "perma"])
def test_disulfide_morph_applies_in_all_modes(mode):
    # disulfide_morph -> apply_bond_morph is the disulfide bonded-term path that
    # replaced the retired alchemy._add_ss_bonded_terms detector.
    t = _Topol()
    apply_bond_morph(t, disulfide_morph(_cys(0), _cys(10), mode=mode, cb="CB2", sg="SG2"))
    assert t.bonds and t.angles and t.dihedrals   # produces the SS bonded set


def test_disulfide_form_A_off_B_on():
    t = _Topol()
    apply_bond_morph(t, disulfide_morph(_cys(0), _cys(10), mode="form", cb="CB2", sg="SG2"))
    # bond: state A force constant zero, state B on
    _, _, a_params, b_params = _sig(t.bonds[0], 2)
    assert a_params[2] == 0.0 and b_params[2] == 100.0


def test_nos_morph_interaction_set():
    cys = {"CA": _At("CA", 0), "CB2": _At("CB2", 1), "SG2": _At("SG2", 2),
           "OD2": _At("OD2", 3)}
    lys = {"CD": _At("CD", 10), "CE2": _At("CE2", 11), "NZ2": _At("NZ2", 12),
           "HZ21": _At("HZ21", 13), "HE21": _At("HE21", 14), "HE22": _At("HE22", 15)}
    m = nos_morph(cys, lys, mode="form")
    kinds = [it.kind for it in m.interactions]
    assert kinds.count("bond") == 1
    assert kinds.count("angle") == 3          # CE-NZ-OD, HZ-NZ-OD, NZ-OD-SG
    assert kinds.count("dihedral") == 6
    # the bridging bond is NZ2-OD2
    bond = next(it for it in m.interactions if it.kind == "bond")
    assert {a.name for a in bond.atoms} == {"NZ2", "OD2"}


def test_nos_morph_optional_atoms_absent():
    """Without HZ21/HE21/HE22 the optional terms are simply dropped."""
    cys = {"CA": _At("CA", 0), "CB2": _At("CB2", 1), "SG2": _At("SG2", 2),
           "OD2": _At("OD2", 3)}
    lys = {"CD": _At("CD", 10), "CE2": _At("CE2", 11), "NZ2": _At("NZ2", 12)}
    m = nos_morph(cys, lys, mode="form")
    kinds = [it.kind for it in m.interactions]
    assert kinds.count("angle") == 2          # no HZ-NZ-OD
    assert kinds.count("dihedral") == 3        # no HZ/HE dihedrals


def test_covalent_edit_container():
    e = CovalentEdit(residues=["C2CD", "C2CD"],
                     bond_morphs=[disulfide_morph(_cys(0), _cys(10), cb="CB2", sg="SG2")])
    assert len(e.bond_morphs) == 1
    assert e.residues == ["C2CD", "C2CD"]
