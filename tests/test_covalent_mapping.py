#!/usr/bin/env python
"""Tests for pmx.covalent.map_residues (graph-MCS core + name prior).

Pure-graph/name logic runs without RDKit-built structures by hand-constructing
small residues; the realistic SER/SEP cases build residues via
TargetResidue.from_rtp and are skipped without RDKit.
"""
import pytest

from pmx.atom import Atom
from pmx.covalent.mapping import _canon_name, map_residues

rdkit = pytest.importorskip("rdkit")  # graph path needs RDKit


class _Res:
    def __init__(self, atoms):
        self.atoms = atoms


def build(spec, bonds):
    """spec: list of (name, element); bonds: list of (name, name)."""
    atoms = [Atom(name=n, symbol=e, id=i + 1) for i, (n, e) in enumerate(spec)]
    by = {a.name: a for a in atoms}
    for a, b in bonds:
        by[a].bonds.append(by[b])
        by[b].bonds.append(by[a])
    return _Res(atoms)


# ── name canonicalization (pure) ──────────────────────────────────────────────

@pytest.mark.parametrize("raw,canon", [
    ("1HB", "HB1"), ("2HB", "HB2"), ("1HG", "HG1"),
    ("HB1", "HB1"), ("HN", "HN"), ("N", "N"), ("3HG2", "HG23"),
])
def test_canon_name(raw, canon):
    assert _canon_name(raw) == canon


# ── graph guarantee: mapping follows topology, not names ──────────────────────

def test_graph_maps_by_position_not_name():
    """A target atom that merely reuses a name at a different graph position
    must NOT capture the WT atom of that name."""
    # WT:     N-CA-CB           ; Target: N-CA-CG-CB  ('CB' now 2 bonds from CA)
    r1 = build([("N", "N"), ("CA", "C"), ("CB", "C")],
               [("N", "CA"), ("CA", "CB")])
    r2 = build([("N", "N"), ("CA", "C"), ("CG", "C"), ("CB", "C")],
               [("N", "CA"), ("CA", "CG"), ("CG", "CB")])
    g = dict(map_residues(r1, r2, method="graph").named_pairs())
    n = dict(map_residues(r1, r2, method="name").named_pairs())
    assert g["CB"] == "CG"     # graph maps by position
    assert n["CB"] == "CB"     # name path is fooled by the label


def test_graph_collision_dummifies_substituents():
    """methyl-CB vs carbonyl-CB: CB<->CB kept via the real CA-CB edge; the
    differing substituents become dummies."""
    r1 = build([("N", "N"), ("CA", "C"), ("CB", "C"),
                ("HB1", "H"), ("HB2", "H"), ("HB3", "H")],
               [("N", "CA"), ("CA", "CB"), ("CB", "HB1"),
                ("CB", "HB2"), ("CB", "HB3")])
    r2 = build([("N", "N"), ("CA", "C"), ("CB", "C"), ("OX", "O")],
               [("N", "CA"), ("CA", "CB"), ("CB", "OX")])
    m = map_residues(r1, r2, method="graph")
    assert dict(m.named_pairs())["CB"] == "CB"
    assert {r1.atoms[i].name for i in m.disappear} == {"HB1", "HB2", "HB3"}
    assert {r2.atoms[j].name for j in m.appear} == {"OX"}


def test_graph_requires_bonds():
    r1 = _Res([Atom(name="N", symbol="N", id=1)])
    r2 = _Res([Atom(name="N", symbol="N", id=1)])
    with pytest.raises(ValueError):
        map_residues(r1, r2, method="graph")


def test_hydrogens_follow_mapped_parent():
    r1 = build([("CA", "C"), ("1HA", "H")], [("CA", "1HA")])
    r2 = build([("CA", "C"), ("HA1", "H")], [("CA", "HA1")])
    m = map_residues(r1, r2, method="graph")
    # heavy CA maps, and its H maps by canonical name across 1HA/HA1
    assert ("CA", "CA") in m.named_pairs()
    assert ("1HA", "HA1") in m.named_pairs()


# ── name-only fallback (unsafe, explicit) ─────────────────────────────────────

def test_name_path_equal_length_no_zip_bug():
    r1 = build([("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O")],
               [("N", "CA"), ("CA", "C"), ("C", "O")])
    r2 = build([("N", "N"), ("HN", "H"), ("CA", "C"), ("HA", "H"),
                ("C", "C"), ("O", "O")],
               [("N", "CA"), ("N", "HN"), ("CA", "C"), ("CA", "HA"), ("C", "O")])
    m = map_residues(r1, r2, method="name")
    assert len(m.pairs1()) == len(m.pairs2())
    assert {r2.atoms[j].name for j in m.appear} == {"HN", "HA"}


# ── realistic SER -> SEP (graph path via from_rtp bonds) ──────────────────────

def _built_residues():
    from pmx.gmx import set_gmxlib
    from pmx.covalent import TargetResidue
    set_gmxlib()
    ser = TargetResidue.from_rtp("SER", "charmm36m-mut").residue
    sep = TargetResidue.from_rtp("SEP", "charmm36m-mut").residue
    return ser, sep


def test_ser_to_sep_graph_maps_hn_and_isolates_phosphate():
    ser, sep = _built_residues()
    m = map_residues(ser, sep)          # auto -> graph (from_rtp carries bonds)
    assert "mcs" in set(m.method.values())
    mapped = dict(m.named_pairs())
    assert mapped["HN"] == "HN" and mapped["OG"] == "OG" and mapped["CB"] == "CB"
    assert {sep.atoms[j].name for j in m.appear} == {"P", "O1P", "O2P", "O3P", "H3T"}
    assert any("HG" in ser.atoms[i].name for i in m.disappear)


def test_proline_backbone_not_decoupled():
    """Proline's ring must not decouple the backbone (backbone anchor overrides
    the ring-only MCS constraint)."""
    from pmx.gmx import set_gmxlib
    from pmx.covalent import TargetResidue
    set_gmxlib()
    ala = TargetResidue.from_rtp("ALA", "amber14sbmut").residue
    pro = TargetResidue.from_rtp("PRO", "amber14sbmut").residue
    m = map_residues(ala, pro)
    mapped = dict(m.named_pairs())
    for bb in ("N", "CA", "C", "O"):
        assert mapped.get(bb) == bb, f"backbone {bb} decoupled"
    # proline ring carbons appear (dummies), they are not mapped to Ala atoms
    assert {pro.atoms[j].name for j in m.appear} >= {"CB", "CG", "CD"}
