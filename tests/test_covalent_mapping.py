#!/usr/bin/env python
"""Tests for pmx.covalent.map_residues (topology/name-first mapper).

The mapper itself needs no RDKit; the realistic SER/SEP integration cases build
their residues via TargetResidue.from_rtp and are skipped without RDKit.
"""
import pytest

from pmx.atom import Atom
from pmx.covalent.mapping import _canon_name, map_residues


class _Res:
    """Minimal residue-like holder for pure-logic tests."""
    def __init__(self, spec):
        # spec: list of (name, element)
        self.atoms = []
        for i, (name, el) in enumerate(spec):
            a = Atom(name=name, symbol=el, id=i + 1)
            self.atoms.append(a)


# ── pure logic (no RDKit) ─────────────────────────────────────────────────────

@pytest.mark.parametrize("raw,canon", [
    ("1HB", "HB1"), ("2HB", "HB2"), ("1HG", "HG1"),
    ("HB1", "HB1"), ("HN", "HN"), ("N", "N"), ("3HG2", "HG23"),
])
def test_canon_name(raw, canon):
    assert _canon_name(raw) == canon


def test_pairs_always_equal_length():
    """Even when res1 lacks atoms present in res2, pairs stay consistent."""
    r1 = _Res([("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O")])  # no H
    r2 = _Res([("N", "N"), ("HN", "H"), ("CA", "C"), ("HA", "H"),
               ("C", "C"), ("O", "O")])
    m = map_residues(r1, r2)
    assert len(m.pairs1()) == len(m.pairs2())
    assert {r1.atoms[i].name for i in m.pairs1()} == {"N", "CA", "C", "O"}
    # res2-only hydrogens are cleanly classified as appearing, not corrupted
    assert {r2.atoms[j].name for j in m.appear} == {"HN", "HA"}
    assert m.disappear == []


def test_name_canonicalization_bridges_pdb_gromacs():
    """1HB (PDB) on one side matches HB1 (GROMACS) on the other."""
    r1 = _Res([("CB", "C"), ("1HB", "H"), ("2HB", "H")])
    r2 = _Res([("CB", "C"), ("HB1", "H"), ("HB2", "H")])
    m = map_residues(r1, r2)
    assert len(m.pairs) == 3
    assert m.disappear == [] and m.appear == []


def test_element_guard_blocks_cross_element_name_collision():
    r1 = _Res([("X1", "C")])
    r2 = _Res([("X1", "O")])
    m = map_residues(r1, r2, element_guard=True)
    assert m.pairs == [] and m.disappear == [0] and m.appear == [0]
    m2 = map_residues(r1, r2, element_guard=False)
    assert len(m2.pairs) == 1


def test_self_map_is_identity():
    r = _Res([("N", "N"), ("CA", "C"), ("CB", "C"), ("OG", "O")])
    m = map_residues(r, _Res([("N", "N"), ("CA", "C"), ("CB", "C"), ("OG", "O")]))
    assert len(m.pairs) == 4
    assert m.disappear == [] and m.appear == []


def test_mcs_fallback_not_implemented():
    r1 = _Res([("N", "N"), ("ZZ", "C")])
    r2 = _Res([("N", "N"), ("QQ", "C")])
    with pytest.raises(NotImplementedError):
        map_residues(r1, r2, mcs_fallback=True)


# ── realistic SER -> SEP (needs RDKit for from_rtp) ───────────────────────────

def _built_residues():
    from pmx.gmx import set_gmxlib
    from pmx.covalent import TargetResidue
    set_gmxlib()
    ser = TargetResidue.from_rtp("SER", "charmm36m-mut").residue
    sep = TargetResidue.from_rtp("SEP", "charmm36m-mut").residue
    return ser, sep


def test_ser_to_sep_maps_hn_by_name():
    pytest.importorskip("rdkit")
    ser, sep = _built_residues()
    m = map_residues(ser, sep)
    mapped_sep = {b for _, b in m.named_pairs()}
    # the whole SER core, hydrogens included, maps -- HN in particular
    assert {"N", "HN", "CA", "HA", "CB", "HB1", "HB2", "OG", "C", "O"} <= mapped_sep
    assert ("HN", "HN") in m.named_pairs()


def test_ser_to_sep_phosphate_appears_hydroxyl_disappears():
    pytest.importorskip("rdkit")
    ser, sep = _built_residues()
    m = map_residues(ser, sep)
    appear = {sep.atoms[j].name for j in m.appear}
    disappear = {ser.atoms[i].name for i in m.disappear}
    assert appear == {"P", "O1P", "O2P", "O3P", "H3T"}
    assert any("HG" in n for n in disappear)
