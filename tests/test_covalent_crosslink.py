#!/usr/bin/env python
"""Tests for pmx.covalent.crosslink -- auto-generated dual-topology cross-links.

Step 1: the sidechain walk that identifies which atoms to duplicate/decouple.
"""
import pytest

pytest.importorskip("rdkit")

from pmx.gmx import set_gmxlib
from pmx.covalent import TargetResidue
from pmx.covalent.crosslink import walk_reacting_atoms

FF = "charmm36m-mut"


@pytest.fixture(scope="module", autouse=True)
def _gmxlib():
    set_gmxlib()


def test_disulfide_walk_matches_c2cd_boundary():
    """CYS/SG depth 2 -> duplicate {CB, SG} + H, boundary CA (as in C2CD)."""
    r = TargetResidue.from_rtp("CYS", FF).residue
    reacting, boundary = walk_reacting_atoms(r, "SG", boundary_depth=2)
    assert boundary == "CA"
    assert set(reacting) == {"CB", "SG", "HB1", "HB2", "HG1"}


def test_nos_lys_walk_matches_k2kn_boundary():
    """LYS/NZ depth 2 -> duplicate {CE, NZ} + H, boundary CD (as in K2KN)."""
    r = TargetResidue.from_rtp("LYS", FF).residue
    reacting, boundary = walk_reacting_atoms(r, "NZ", boundary_depth=2)
    assert boundary == "CD"
    assert {"CE", "NZ"} <= set(reacting)
    assert "CD" not in reacting            # boundary stays common


def test_walk_depth_one_stops_earlier():
    r = TargetResidue.from_rtp("CYS", FF).residue
    reacting, boundary = walk_reacting_atoms(r, "SG", boundary_depth=1)
    assert boundary == "CB"
    assert set(reacting) == {"SG", "HG1"}


def test_walk_unknown_bond_atom_raises():
    r = TargetResidue.from_rtp("CYS", FF).residue
    with pytest.raises(KeyError):
        walk_reacting_atoms(r, "ZZ")


def test_build_split_residue_reproduces_c2cd_structure():
    """Auto-generated CYS+CYS2 split-atom hybrid matches C2CD structure/params."""
    from pmx.covalent.crosslink import build_split_residue
    cys = TargetResidue.from_rtp("CYS", FF).residue
    cys2 = TargetResidue.from_rtp("CYS2", FF).residue
    atoms, bonds = build_split_residue(cys, cys2, "SG", boundary_depth=2)
    by = {a.name: a for a in atoms}

    assert len(atoms) == 15                         # same count as C2CD

    # free branch (real A / dummy B) with WT free-thiol params
    assert by["SG1"].atomtype == "S" and by["SG1"].atomtypeB == "DUM_S"
    assert by["CB1"].atomtype == "CT2" and by["CB1"].atomtypeB == "DUM_CT2"
    # bound branch (dummy A / real B) with disulfide params
    assert by["SG2"].atomtype == "DUM_SM" and by["SG2"].atomtypeB == "SM"
    assert by["CB2"].atomtype == "DUM_CT2" and by["CB2"].atomtypeB == "CT2"
    # thiol H lives on the free branch only and disappears
    assert "HG11" in by and by["HG11"].atomtypeB.startswith("DUM")
    assert not any(n.startswith("HG") and n != "HG11" for n in by)

    # boundary CA (common) reconnects to BOTH branches
    assert ("CA", "CB1") in bonds and ("CA", "CB2") in bonds
    # backbone stays common (single, morphing atoms)
    for bb in ("N", "CA", "C", "O"):
        assert bb in by and not by[bb].atomtype.startswith("DUM")
