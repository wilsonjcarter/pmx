#!/usr/bin/env python
"""Tests for pmx.covalent.TargetResidue (RTP-name target front end).

Skipped automatically when RDKit is not installed.
"""
import math
import os
from functools import lru_cache

import pytest

rdkit = pytest.importorskip("rdkit")

from pmx.gmx import set_gmxlib
from pmx.model import Model
from pmx.ffparser import RTPParser
from pmx.utils import get_ff_path
from pmx.covalent import TargetResidue

FF = "charmm36m-mut"


@pytest.fixture(scope="module", autouse=True)
def _gmxlib():
    # point get_ff_path at the pmx built-in mutff force fields
    set_gmxlib()


def _dist(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


@lru_cache(maxsize=4)
def _rtp_oracle(ff):
    ffpath = get_ff_path(ff)
    f = os.path.join(ffpath, "merged.rtp")
    if not os.path.isfile(f):
        f = os.path.join(ffpath, "aminoacids.rtp")
    return RTPParser(f)


def _rtp_reference(resname, ff):
    entry = _rtp_oracle(ff)[resname]
    names = [a[0] for a in entry["atoms"]]
    bonds = [(a1, a2) for a1, a2 in entry["bonds"]
             if not a1.startswith(("+", "-")) and not a2.startswith(("+", "-"))]
    return names, bonds


# Representative chemistries: aliphatic, hydroxyl, aromatic, fused-ring,
# cationic (Lys/Arg/His+), anionic (Asp/Glu), ring-on-backbone (Pro),
# and phospho-PTMs (the headline covalent targets).
_RESIDUES = ["SER", "THR", "TYR", "TRP", "LYS", "ARG", "HSP",
             "ASP", "GLU", "PRO", "GLY", "SEP", "TPO", "PTR"]


@pytest.fixture(scope="module")
def built(_gmxlib):
    """Build each target residue once (embedding + RTP parse are not cheap)."""
    out = {}
    for r in _RESIDUES:
        tr = TargetResidue.from_rtp(r, FF)
        ref_names, ref_bonds = _rtp_reference(r, FF)
        out[r] = (tr, ref_names, ref_bonds)
    return out


@pytest.mark.parametrize("resname", _RESIDUES)
def test_from_rtp_names_match(built, resname):
    """Built residue has exactly the RTP atoms, in order."""
    tr, ref_names, _ = built[resname]
    assert [a.name for a in tr.residue.atoms] == ref_names


@pytest.mark.parametrize("resname", _RESIDUES)
def test_from_rtp_elements_assigned(built, resname):
    """Every atom gets a non-empty element symbol; H/C/N/O present as expected."""
    tr, _, _ = built[resname]
    symbols = [a.symbol for a in tr.residue.atoms]
    assert all(symbols), "some atom has empty element symbol"
    # every standard residue has carbon and hydrogen
    assert "C" in symbols and "H" in symbols


@pytest.mark.parametrize("resname", _RESIDUES)
def test_from_rtp_geometry_sane(built, resname):
    """Embedded geometry is non-degenerate with chemical bond lengths."""
    tr, _, ref_bonds = built[resname]
    by_name = {a.name: a for a in tr.residue.atoms}
    for n1, n2 in ref_bonds:
        d = _dist(by_name[n1].x, by_name[n2].x)
        assert 0.85 <= d <= 2.0, f"{resname} {n1}-{n2} = {d:.2f} A out of range"


def test_from_rtp_pdb_roundtrip(tmpdir):
    """Target writes a PDB that reloads with the same atom count."""
    tr = TargetResidue.from_rtp("SEP", FF)
    out = str(tmpdir.join("sep_target.pdb"))
    tr.write(out)
    assert len(Model(out).atoms) == len(tr.residue.atoms)


def test_from_rtp_deterministic():
    """Fixed seed gives identical coordinates across builds."""
    a = TargetResidue.from_rtp("TYR", FF, seed=123)
    b = TargetResidue.from_rtp("TYR", FF, seed=123)
    for x, y in zip(a.residue.atoms, b.residue.atoms):
        assert x.x == pytest.approx(y.x, abs=1e-6)


def test_from_rtp_missing_residue():
    """A residue absent from the RTP raises a clear KeyError."""
    with pytest.raises(KeyError):
        TargetResidue.from_rtp("ZZZ_NOT_A_RESIDUE", FF)
