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


def _two_residue_model(resname, resids=(1, 2), gap=6.0):
    """Build a minimal 2-residue Model from a force-field RTP residue.

    Writes the residue's atoms twice (translated by *gap*) to a temporary PDB
    and reads it back -- enough for the model-level cross-link placement path
    without needing gmx/pdb2gmx.
    """
    import os
    import tempfile
    import numpy as np
    from pmx.model import Model

    res = TargetResidue.from_rtp(resname, FF).residue

    def line(serial, name, resid, x):
        nm = name if len(name) == 4 else " " + name
        return (f"ATOM  {serial:>5} {nm:<4} {resname:>3} A{resid:>4}    "
                f"{x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}  1.00  0.00\n")

    lines, serial = [], 1
    for i, resid in enumerate(resids):
        shift = np.array([gap * i, 0.0, 0.0])
        for a in res.atoms:
            lines.append(line(serial, a.name, resid, np.array(a.x) + shift))
            serial += 1
    lines.append("TER\nEND\n")
    fn = tempfile.mktemp(suffix=".pdb")
    with open(fn, "w") as fh:
        fh.writelines(lines)
    try:
        return Model(fn, bPDBTER=True, rename_atoms=True)
    finally:
        os.remove(fn)


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


def test_emit_crosslink_rtp_builds_c2cd_entry():
    """emit_crosslink_rtp produces a valid split-atom RTP entry (for pdb2gmx)."""
    from pmx.covalent.crosslink import decouple, emit_crosslink_rtp
    cys = TargetResidue.from_rtp("CYS", FF).residue
    cys2 = TargetResidue.from_rtp("CYS2", FF).residue
    _amap, atoms, bonds = decouple(cys, cys2, "SG")
    rtp = emit_crosslink_rtp(atoms, bonds, "C2CD", "CYS", "CYS2", FF)
    assert "[ C2CD ]" in rtp
    for sec in ("[ atoms ]", "[ bonds ]", "[ impropers ]"):
        assert sec in rtp
    # split branches present in the atom block
    for n in ("CB1", "CB2", "SG1", "SG2"):
        assert n in rtp


def test_split_branches_colocated_on_wt():
    """Each split atom (free X1 / bound X2) sits exactly on its WT counterpart.

    The bound branch represents the same physical atom as the free branch, so
    both must reuse the WT coordinate -- otherwise the classic split placement
    (which co-locates CB1/CB2, SG1/SG2) is not reproduced.
    """
    import numpy as np
    from pmx.covalent.crosslink import build_split_residue
    cys = TargetResidue.from_rtp("CYS", FF).residue
    cys2 = TargetResidue.from_rtp("CYS2", FF).residue
    wt = {a.name: np.array(a.x) for a in cys.atoms}
    atoms, _ = build_split_residue(cys, cys2, "SG")
    by = {a.name: np.array(a.x) for a in atoms}
    for base in ("CB", "SG"):
        assert np.allclose(by[base + "1"], by[base + "2"])   # branches coincident
        assert np.allclose(by[base + "1"], wt[base])         # and on the WT atom


def test_decouple_mutate_rejects_bad_react_atom():
    """A reacting atom absent from the residue is a clear ValueError."""
    from pmx.covalent import decouple_mutate
    m = _two_residue_model("CYS")
    with pytest.raises(ValueError, match="no atom 'ZZ'"):
        decouple_mutate(m, 1, "CYS2", FF, react_atom="ZZ",
                        inplace=True, write_files=False)


def test_decouple_mutate_rejects_missing_target():
    """A cross-link target not in the force field points at endpoint construction."""
    from pmx.covalent import decouple_mutate
    m = _two_residue_model("CYS")
    with pytest.raises(ValueError, match="not in force field"):
        decouple_mutate(m, 1, "NOPE", FF, react_atom="SG",
                        inplace=True, write_files=False)


def test_decouple_mutate_rejects_atom_missing_from_target():
    """The reacting atom must have a bound-state counterpart to split into."""
    from pmx.covalent import decouple_mutate
    m = _two_residue_model("CYS")
    with pytest.raises(ValueError, match="absent from the bound target"):
        decouple_mutate(m, 1, "ALA", FF, react_atom="SG",
                        inplace=True, write_files=False)


def test_split_bonds_names_bound_only_bridging_atom():
    """A target-only bridging atom (e.g. a NOS/sulfenyl OD) is a branch-2 atom.

    Its bonds must carry the branch-2 suffix, and a bond between two bound-only
    atoms (OD-HD) must not be dropped -- otherwise the emitted connectivity
    references a nonexistent 'OD' or loses the O-H bond.
    """
    from pmx.covalent.crosslink import decouple
    cys = TargetResidue.from_rtp("CYS", FF).residue
    cso = TargetResidue.from_rtp("CSO", FF).residue   # Cys-sulfenic: SG-OD-HD
    _amap, atoms, bonds = decouple(cys, cso, "SG")
    names = {a.name for a in atoms}
    assert "OD2" in names and "HD2" in names          # bound-only atoms suffixed
    assert not any(b[0] not in names or b[1] not in names for b in bonds)
    assert any(set(b) == {"OD2", "SG2"} for b in bonds)   # sulfenyl bond named
    assert any(set(b) == {"OD2", "HD2"} for b in bonds)   # O-H bond not dropped


# -- cross-link hybrid naming (regression: LYS must be K2KN, not L2KN) ----------

def test_nos_lys_hybrid_named_k2kn():
    """Lys NOS endpoint must use the one-letter code K, giving K2KN.

    Regression for the resnA[:1] bug: the resname 'LYS' starts with 'L', so the
    old code produced 'L2KN' -- unregistered in residuetypes.dat and unrouted by
    get_mtp_file (both expect K2KN)."""
    from pmx.covalent import decouple_mutate
    m = _two_residue_model("LYS")
    name = decouple_mutate(m, 1, "KN0", FF, react_atom="NZ", write_files=False)
    assert name == "K2KN", name


def test_nos_cys_hybrid_named_c2cn():
    from pmx.covalent import decouple_mutate
    m = _two_residue_model("CYS")
    name = decouple_mutate(m, 1, "CN0", FF, react_atom="SG", write_files=False)
    assert name == "C2CN", name
