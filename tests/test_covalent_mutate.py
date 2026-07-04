#!/usr/bin/env python
"""Tests for pmx.covalent.build_hybrid_residue / mutate (graph-mapper emit).

Verifies the emitted hybrid MTP maps shared atoms (notably the polar HN) by
identity rather than dropping them to dummies -- the defect of the domapping +
polar-H-patch path. Needs RDKit (for from_rtp) and the pmx force fields.
"""
import pytest

pytest.importorskip("rdkit")

from pmx.gmx import set_gmxlib
from pmx.covalent import TargetResidue, build_hybrid_residue

FF = "charmm36m-mut"


def _morphes(mtp_text):
    """Parse the [ morphes ] section into (a_name, a_type, b_name, b_type)."""
    out, inseg = [], False
    for line in mtp_text.splitlines():
        s = line.strip()
        if s.startswith("[ morphes"):
            inseg = True
            continue
        if inseg and s.startswith("["):
            break
        if inseg and "->" in s:
            left, right = s.split("->")
            la = left.split()
            ra = right.split()
            if len(la) >= 2 and len(ra) >= 2:
                out.append((la[0], la[1], ra[0], ra[1]))
    return out


@pytest.fixture(scope="module")
def hybrid():
    set_gmxlib()
    ser = TargetResidue.from_rtp("SER", FF).residue
    sep = TargetResidue.from_rtp("SEP", FF).residue
    rtp_text, mtp_text, name = build_hybrid_residue(
        ser, sep, FF, "SEP", write_ff_types=False)
    return rtp_text, mtp_text, name, _morphes(mtp_text)


def test_hybrid_name(hybrid):
    assert hybrid[2] == "S2SE"


def test_hn_is_common_not_dummy(hybrid):
    """HN must map HN->HN with no DUM_ type on either side (the HN-bug fix)."""
    morphes = hybrid[3]
    hn = [m for m in morphes if m[0] == "HN"]
    assert hn, "HN missing from morphes"
    a_name, a_type, b_name, b_type = hn[0]
    assert b_name == "HN"
    assert not a_type.startswith("DUM") and not b_type.startswith("DUM")


def test_backbone_and_sidechain_core_common(hybrid):
    """The whole SER core maps by identity (no DUM_ on either state)."""
    morphes = {m[0]: m for m in hybrid[3]}
    for name in ("N", "HN", "CA", "HA", "CB", "OG", "C", "O"):
        assert name in morphes, f"{name} not in morphes"
        _, a_type, _, b_type = morphes[name]
        assert not a_type.startswith("DUM"), f"{name} is a state-A dummy"


def test_hydroxyl_h_disappears(hybrid):
    """The SER hydroxyl H is real in A and dummy in B (.gone)."""
    morphes = hybrid[3]
    gone = [m for m in morphes if "gone" in m[2] or m[3].startswith("DUM")]
    assert any("HG" in m[0] for m in gone), "hydroxyl H should disappear"


def test_mutate_from_structure_matches_rtp_path(tmp_path):
    """structure= (geometry from a PDB) yields the same hybrid as the default
    RTP path -- only coordinates differ."""
    from pmx.gmx import set_gmxlib
    from pmx.model import Model
    from pmx.covalent import TargetResidue, mutate
    set_gmxlib()
    ser = str(tmp_path / "ser.pdb")
    sep = str(tmp_path / "sep.pdb")
    TargetResidue.from_rtp("SER", FF).write(ser)      # protonated WT
    TargetResidue.from_rtp("SEP", FF).write(sep)      # target geometry

    m_def = Model(ser)
    n_def = mutate(m_def, 1, "SEP", FF, inplace=True, write_files=False)
    m_str = Model(ser)
    n_str = mutate(m_str, 1, "SEP", FF, structure=sep,
                   inplace=True, write_files=False)

    assert n_def == n_str == "S2SE"
    assert sorted(a.name for a in m_def.residues[0].atoms) == \
        sorted(a.name for a in m_str.residues[0].atoms)
    # HN common (backbone) and phosphate present as dummies
    assert any(a.name == "HN" for a in m_str.residues[0].atoms)
    assert any(a.name == "DP" for a in m_str.residues[0].atoms)


def test_itp_target_requires_structure():
    from pmx.covalent.mutate import _build_target_residue
    with pytest.raises(ValueError):
        _build_target_residue("SEP", FF, {}, None, None,
                              structure=None, itp="frag.itp")


def _write_capped_sep_fragment(tmp_path):
    """Build a CH3-capped SEP sidechain fragment (ITP + PDB) from from_rtp SEP."""
    import numpy as np
    from pmx.covalent import TargetResidue
    keep = ["CA", "HA", "CB", "HB1", "HB2", "OG", "P", "O1P", "O2P", "O3P", "H3T"]
    mass = {"C": 12.011, "H": 1.008, "O": 15.999, "P": 30.974}
    sep = TargetResidue.from_rtp("SEP", FF).residue
    by = {a.name: a for a in sep.atoms}
    ca, cb = np.array(by["CA"].x), np.array(by["CB"].x)
    d = (ca - cb) / np.linalg.norm(ca - cb)
    cy = ca + d * 1.5
    hy = [cy + np.array(v) for v in ([.6, .6, 0], [-.6, .6, 0], [0, -.6, .6])]
    frag = [(n, by[n].atomtype, by[n].q, (by[n].symbol or n[0]), by[n].x)
            for n in keep]
    frag += [("CY", "CT3", -0.27, "C", list(cy))]
    frag += [(f"HY{i}", "HA3", 0.09, "H", list(h)) for i, h in enumerate(hy, 1)]
    itp, pdb = str(tmp_path / "f.itp"), str(tmp_path / "f.pdb")
    nr = {n: i for i, (n, *_ ) in enumerate(frag, 1)}
    with open(itp, "w") as f:
        f.write("[ moleculetype ]\nSEP 3\n\n[ atoms ]\n")
        for i, (n, t, q, el, x) in enumerate(frag, 1):
            f.write(f"{i} {t} 1 SEP {n} {i} {q:.4f} {mass[el[0]]:.4f}\n")
        f.write("\n[ bonds ]\n")
        for a, b in [("CA", "HA"), ("CA", "CB"), ("CB", "HB1"), ("CB", "HB2"),
                     ("CB", "OG"), ("OG", "P"), ("P", "O1P"), ("P", "O2P"),
                     ("P", "O3P"), ("O3P", "H3T"), ("CA", "CY"), ("CY", "HY1"),
                     ("CY", "HY2"), ("CY", "HY3")]:
            f.write(f"{nr[a]} {nr[b]} 1\n")
    with open(pdb, "w") as f:
        for i, (n, t, q, el, x) in enumerate(frag, 1):
            f.write(f"ATOM  {i:5} {n:<4} SEP A   1    "
                    f"{x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}  1.00  0.00          "
                    f"{el[0]:>2}\n")
        f.write("END\n")
    return itp, pdb


def test_mutate_via_itp_fragment(tmp_path):
    """End-to-end itp path: capped fragment -> hybrid with phosphate dummies."""
    from pmx.gmx import set_gmxlib
    from pmx.model import Model
    from pmx.covalent import TargetResidue, mutate
    set_gmxlib()
    itp, pdb = _write_capped_sep_fragment(tmp_path)
    TargetResidue.from_rtp("SER", FF).write(str(tmp_path / "wt.pdb"))
    m = Model(str(tmp_path / "wt.pdb"))
    name = mutate(m, 1, "SEP", FF, itp=itp, structure=pdb,
                  inplace=True, write_files=False)
    names = [a.name for a in m.residues[0].atoms]
    assert name == "S2SE"
    assert "HN" in names                      # WT backbone H mapped through
    assert any(a.name == "DP" for a in m.residues[0].atoms)  # phosphate appears


def test_protonation_warning_on_hydrogenless_wt(tmp_path):
    from pmx.gmx import set_gmxlib
    from pmx.model import Model
    from pmx.covalent import TargetResidue, mutate
    set_gmxlib()
    tr = TargetResidue.from_rtp("SER", FF)
    tr.residue.atoms = [a for a in tr.residue.atoms if a.symbol != "H"]
    tr.model.atoms = tr.residue.atoms
    p = str(tmp_path / "ser_heavy.pdb")
    tr.model.write(p)
    m = Model(p)
    with pytest.warns(UserWarning, match="unprotonated"):
        mutate(m, 1, "SEP", FF, inplace=True, write_files=False)


def test_cys_thiol_to_disulfide_type_morph_maps():
    """The graph core maps the CYS -> disulfide-Cys charge/type change: SG
    retypes (S -> SM, free thiol -> disulfide) and the thiol hydrogen disappears.

    NOTE: this is ONLY the per-residue electrostatic/type morph.  A physically
    valid alchemical disulfide additionally needs the SG-SG bond's topological
    1-4 exclusions handled (they would otherwise corrupt the free-thiol state) --
    via the vsite + [ pairs ] scheme documented in
    docs/developers/covalent_framework.md.  A plain 1:1 residue map does NOT
    provide that, so this test asserts the mapping only, not a valid disulfide.
    """
    from pmx.gmx import set_gmxlib
    from pmx.covalent import TargetResidue, build_hybrid_residue, map_residues
    set_gmxlib()
    cys = TargetResidue.from_rtp("CYS", FF).residue
    cys2 = TargetResidue.from_rtp("CYS2", FF).residue        # disulfide cystine
    m = map_residues(cys, cys2)
    assert dict(m.named_pairs())["SG"] == "SG"               # SG common (retypes)
    assert any("HG" in cys.atoms[i].name for i in m.disappear)   # thiol H goes
    _, mtp, _ = build_hybrid_residue(cys, cys2, FF, "CYS2", write_ff_types=False)
    sg = [l.strip() for l in mtp.splitlines()
          if l.strip().startswith("SG ") and "->" in l]
    assert sg and "SM" in sg[0] and "DUM" not in sg[0]


def test_split_rtp_header():
    from pmx.covalent.mutate import _split_rtp_header
    rtp = "[ bondedtypes ]\nHDR\n\n[ A2SE ] ; x\n[ atoms ]\n  N N\n"
    h, e = _split_rtp_header(rtp, "A2SE")
    assert h == "[ bondedtypes ]\nHDR\n\n"
    assert e.startswith("[ A2SE ]")


def test_upsert_entry_create_append_replace(tmp_path):
    from pmx.covalent.mutate import _upsert_entry
    p = str(tmp_path / "nsaa.rtp")
    # create with header
    _upsert_entry(p, "A2SE", "[ A2SE ] ; a\n[ atoms ]\n  N N\n",
                  header="[ bondedtypes ]\nHDR\n\n")
    txt = open(p).read()
    assert txt.startswith("[ bondedtypes ]")
    assert txt.count("[ A2SE ]") == 1
    # append a different hybrid -- header not duplicated
    _upsert_entry(p, "B2SE", "[ B2SE ] ; b\n[ atoms ]\n  N N\n")
    txt = open(p).read()
    assert "[ A2SE ]" in txt and "[ B2SE ]" in txt
    assert txt.count("[ bondedtypes ]") == 1
    # re-mutate A2SE -> its entry is replaced, not duplicated; B2SE untouched
    _upsert_entry(p, "A2SE", "[ A2SE ] ; a2\n[ atoms ]\n  N NEW\n")
    txt = open(p).read()
    assert txt.count("[ A2SE ]") == 1
    assert "NEW" in txt and "; a2" in txt
    assert "[ B2SE ]" in txt


def test_phosphate_appears_as_bstate(hybrid):
    """Phosphate atoms are dummies in A (DUM_) and real in B."""
    morphes = {m[2]: m for m in hybrid[3]}  # keyed by B-state name
    for bname in ("P", "O1P", "O2P", "O3P"):
        assert bname in morphes, f"{bname} not appearing in B state"
        a_name, a_type, _, b_type = morphes[bname]
        assert a_type.startswith("DUM"), f"{bname} should be an A-state dummy"
        assert not b_type.startswith("DUM")
