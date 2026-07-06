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


def test_build_crosslink_mtp_emits_c2cd_morphes():
    """The split-atom hybrid emits a valid MTP with the C2CD morph structure."""
    from pmx.covalent.crosslink import build_crosslink_mtp
    cys = TargetResidue.from_rtp("CYS", FF).residue
    cys2 = TargetResidue.from_rtp("CYS2", FF).residue
    mtp = build_crosslink_mtp(cys, cys2, "SG", "C2CD")
    for section in ("[ morphes ]", "[ atoms ]", "[ coords ]"):
        assert section in mtp
    morphes = {l.split()[0]: l for l in mtp.splitlines()
               if l.strip() and l.split()[0] in
               ("SG1", "SG2", "CB1", "CB2")}
    # free branch disappears (A-real type -> DUM), bound branch appears
    assert "S" in morphes["SG1"] and "DUM_S" in morphes["SG1"]
    assert "DUM_SM" in morphes["SG2"] and "SM" in morphes["SG2"]
    assert "CT2" in morphes["CB1"] and "DUM_CT2" in morphes["CB1"]
    # backbone impropers emitted (as in C2CD)
    assert "[ impropers ]" in mtp
    assert "N      -C     CA     HN" in mtp and "C      CA     +N     O" in mtp


class _BP:
    def get_bond_param(self, *t):     return (1, 2.05, 300.0)
    def get_angle_param(self, *t):    return (5, 110.0, 50.0, 0.0, 0.0)
    def get_dihedral_param(self, *t): return [(9, 180.0, 2.0, 2)]


class _Topol:
    def __init__(self):
        self.bonds, self.angles, self.dihedrals = [], [], []
        self.BondedParams = _BP()


def test_disulfide_edit_composes_decouple_and_bond():
    """A full disulfide = decouple x2 + bond, on the bound branch."""
    from pmx.covalent.crosslink import disulfide_edit
    from pmx.covalent.edit import apply_bond_morph
    cys = lambda: TargetResidue.from_rtp("CYS", FF).residue      # noqa: E731
    cys2 = lambda: TargetResidue.from_rtp("CYS2", FF).residue    # noqa: E731
    edit, s1, s2 = disulfide_edit(cys(), cys(), cys2(), cys2())

    assert len(edit.bond_morphs) == 1
    morph = edit.bond_morphs[0]
    bond = next(it for it in morph.interactions if it.kind == "bond")
    # the SG-SG bond sits on the bound-branch SG2 (dummy in A)
    assert all(a.name == "SG2" and a.atomtype == "DUM_SM" for a in bond.atoms)

    t = _Topol()
    apply_bond_morph(t, morph)               # form mode
    assert t.bonds[0][3][2] == 0.0           # state A off
    assert t.bonds[0][4][2] == 300.0         # state B on
    assert len(t.bonds) + len(t.angles) + len(t.dihedrals) == 6


def test_composable_covalent_edit_api():
    """Compose a disulfide from the first-class operations: decouple x2 + bond."""
    from pmx.covalent import CovalentEdit
    from pmx.covalent.edit import disulfide_morph
    cys = lambda: TargetResidue.from_rtp("CYS", FF).residue      # noqa: E731
    cys2 = lambda: TargetResidue.from_rtp("CYS2", FF).residue    # noqa: E731

    edit = CovalentEdit()
    a1 = edit.decouple(cys(), cys2(), "SG")     # returns name->Atom
    a2 = edit.decouple(cys(), cys2(), "SG")
    edit.bond(disulfide_morph(a1, a2, cb="CB2", sg="SG2"))

    assert len(edit.decoupled) == 2 and len(edit.bond_morphs) == 1
    assert edit.residues == ["CYS", "CYS"]
    # decouple exposes the bound-branch SG2 (dummy in A) for the bond
    assert a1["SG2"].atomtype == "DUM_SM"


def test_covalent_edit_unified_emit_single_and_crosslink():
    """A single-residue mutation and a cross-link emit through one CovalentEdit."""
    from pmx.covalent import CovalentEdit
    # single-residue mutation via the same API
    e1 = CovalentEdit(ff=FF).mutate(
        TargetResidue.from_rtp("SER", FF).residue,
        TargetResidue.from_rtp("SEP", FF).residue)
    assert "[ morphes ]" in e1.emit() and e1.hybrids[0][0] == "S2SE"

    # cross-link (decoupled) via the same API + same .emit()
    e2 = CovalentEdit()
    e2.decouple(TargetResidue.from_rtp("CYS", FF).residue,
                TargetResidue.from_rtp("CYS2", FF).residue, "SG")
    m2 = e2.emit()
    assert "SG1" in m2 and "SG2" in m2


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


def test_crosslink_places_split_residue_in_model():
    """Model-level crosslink: both cysteines become the split hybrid in place."""
    import numpy as np
    from pmx.covalent import crosslink
    m = _two_residue_model("CYS")
    # write_files=False -> no force-field writes (side-effect free test)
    name = crosslink(m, 1, 2, FF, target1="CYS2", inplace=True,
                     write_files=False)
    assert name == "C2CY"
    for resid in (1, 2):
        r = m.fetch_residue(idx=resid)
        assert r.resname == name
        by = {a.name: np.array(a.x) for a in r.atoms}
        for n in ("CB1", "CB2", "SG1", "SG2", "HG11"):
            assert n in by, f"resid {resid} missing {n}"
        # split branches co-located; state-A masses assigned (not zero)
        assert np.allclose(by["SG1"], by["SG2"])
        assert np.allclose(by["CB1"], by["CB2"])
        assert all(a.m > 0 for a in r.atoms)


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


def test_crosslink_hetero_emits_two_entries(tmp_path):
    """A hetero cross-link (NOS-shaped: Cys vs Lys) places two distinct hybrids,
    with the forming bond declared once via an explicit link atom.

    Uses CSO/LYN as bound-target proxies to exercise the machinery (two entries,
    link_atom != reacting atom, bond on residue 1); the chemically-correct NOS
    target is a separate force-field data question.
    """
    from pmx.covalent import crosslink
    # build a CYS + LYS model on the fly
    import os
    import numpy as np
    from pmx.model import Model
    cys = TargetResidue.from_rtp("CYS", FF).residue
    lys = TargetResidue.from_rtp("LYS", FF).residue

    def block(res, resname, resid, shift, start):
        out = []
        for i, a in enumerate(res.atoms):
            nm = a.name if len(a.name) == 4 else " " + a.name
            x = np.array(a.x) + shift
            out.append(f"ATOM  {start + i:>5} {nm:<4} {resname:>3} A{resid:>4}    "
                       f"{x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}  1.00  0.00\n")
        return out

    lines = block(cys, "CYS", 1, np.zeros(3), 1)
    lines += block(lys, "LYS", 2, np.array([8.0, 0, 0]), len(cys.atoms) + 1)
    lines.append("TER\nEND\n")
    fn = str(tmp_path / "cyslys.pdb")
    with open(fn, "w") as fh:
        fh.writelines(lines)
    m = Model(fn, bPDBTER=True, rename_atoms=True)

    crosslink(m, 1, 2, FF, target1="CSO", target2="LYN",
              bond_atom1="SG", bond_atom2="NZ",
              link_atom1="OD", link_atom2="NZ",
              inplace=True, write_files=False)
    n1 = m.fetch_residue(idx=1).resname
    n2 = m.fetch_residue(idx=2).resname
    assert n1 != n2                                    # hetero: distinct hybrids
    by1 = {a.name for a in m.fetch_residue(idx=1).atoms}
    by2 = {a.name for a in m.fetch_residue(idx=2).atoms}
    assert "OD2" in by1                                # Cys bridging atom placed
    assert "NZ2" in by2                                # Lys reacting atom placed
