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
