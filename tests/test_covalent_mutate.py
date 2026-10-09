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


def _atom_charges(mtp_text):
    """Parse the [ atoms ] section into (chargeA, chargeB) per atom."""
    out, inseg = [], False
    for line in mtp_text.splitlines():
        s = line.split(";")[0].strip()
        if s.startswith("[ atoms"):
            inseg = True
            continue
        if inseg and s.startswith("["):
            break
        f = s.split()
        # name typeA chargeA cgnr massA typeB chargeB massB
        if inseg and len(f) >= 7:
            try:
                out.append((float(f[2]), float(f[6])))
            except ValueError:
                pass
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
    assert hybrid[2] == "S2SP"  # SEP first+last (was truncated S2SE)


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


def _tripeptide_model(resname):
    """A 3-residue chain of *resname* so residue 2 is internal (not a terminus)."""
    from copy import deepcopy
    from pmx.model import Model
    base = TargetResidue.from_rtp(resname, FF).residue
    atoms, aid = [], 1
    for ri in (1, 2, 3):
        for a in base.atoms:
            b = deepcopy(a)
            b.id = aid
            b.resnr = ri
            b.chain_id = "A"
            b.x = [a.x[0] + ri * 5.0, a.x[1], a.x[2]]
            atoms.append(b)
            aid += 1
    return Model(atoms=atoms, renumber_residues=False)


def test_mutate_from_structure_matches_rtp_path(tmp_path):
    """structure= (geometry from a PDB) yields the same hybrid as the default
    RTP path -- only coordinates differ.  Uses an internal residue (residue 2 of
    a tripeptide) so no terminal patch is involved."""
    from pmx.gmx import set_gmxlib
    from pmx.covalent import TargetResidue, mutate
    set_gmxlib()
    sep = str(tmp_path / "sep.pdb")
    TargetResidue.from_rtp("SEP", FF).write(sep)      # target geometry

    m_def = _tripeptide_model("SER")
    n_def = mutate(m_def, 2, "SEP", FF, inplace=True, write_files=False)
    m_str = _tripeptide_model("SER")
    n_str = mutate(m_str, 2, "SEP", FF, structure=sep,
                   inplace=True, write_files=False)

    assert n_def == n_str == "S2SP"
    assert sorted(a.name for a in m_def.residues[1].atoms) == \
        sorted(a.name for a in m_str.residues[1].atoms)
    # HN common (backbone) and phosphate present as dummies
    assert any(a.name == "HN" for a in m_str.residues[1].atoms)
    assert any(a.name == "DP" for a in m_str.residues[1].atoms)


def test_terminal_mutation_routes_generative_and_conserves_charge(tmp_path):
    """A residue at a terminus is mutated through the generative engine with the
    terminus patched in both states -- a charge-neutral swap stays neutral."""
    from pmx.gmx import set_gmxlib
    from pmx.covalent import build_hybrid_residue
    set_gmxlib()
    # C-terminal ALA -> SER: both neutral, so state-A and state-B residue
    # charges must match when both are built in terminal form.
    r1 = TargetResidue.from_rtp("ALA", FF, term="C").residue
    r2 = TargetResidue.from_rtp("SER", FF, term="C").residue
    _, mtp, hn = build_hybrid_residue(r1, r2, FF, "SER", term="C",
                                      write_ff_types=False)
    qA = sum(m[0] for m in _atom_charges(mtp))
    qB = sum(m[1] for m in _atom_charges(mtp))
    assert round(qA) == round(qB)          # neutral -> neutral at the terminus
    # the C-terminal carboxylate (OT2) is present, not an internal single O
    names = [m[0] for m in _morphes(mtp)]
    assert "OT2" in names or any(n.startswith("OT") for n in names)


def test_terminal_neighbour_morph_cterm():
    # A terminal deletion's neighbour is an internal->terminal residue morph:
    # build_hybrid_residue(internal, terminal) via the same engine as any
    # mutation, producing the C-terminus cap (C->CC, O->OT1, appearing OT2).
    set_gmxlib()
    r1 = TargetResidue.from_rtp("VAL", FF).residue
    r2 = TargetResidue.from_rtp("VAL", FF, term="C").residue
    _, mtp, hn = build_hybrid_residue(r1, r2, FF, "VAL", term=None,
                                      term_b="C", hybrid_name="VdeC",
                                      write_ff_types=False)
    assert hn == "VdeC"
    morphs = {m[0]: m for m in _morphes(mtp)}   # keyed by state-A atom name
    # carbonyl C retypes C -> CC; O is renamed OT1 and retyped O -> OC
    assert morphs["C"][3] == "CC"
    assert morphs["O"][2] == "OT1" and morphs["O"][3] == "OC"
    # the second carboxyl oxygen appears from a dummy in state B
    appearing = [m for m in _morphes(mtp) if m[1].startswith("DUM")]
    assert any(m[3] == "OC" for m in appearing)


def test_terminal_neighbour_morph_nterm_adds_two_protons():
    set_gmxlib()
    r1 = TargetResidue.from_rtp("ALA", FF).residue
    r2 = TargetResidue.from_rtp("ALA", FF, term="N").residue
    _, mtp, hn = build_hybrid_residue(r1, r2, FF, "ALA", term=None,
                                      term_b="N", hybrid_name="AdeN",
                                      write_ff_types=False)
    assert hn == "AdeN"
    morphs = {m[0]: m for m in _morphes(mtp)}
    assert morphs["N"][3] == "NH3"          # amide N -> charged N-terminus
    # two amide protons appear from dummies (H2, H3)
    appearing = [m for m in _morphes(mtp) if m[1].startswith("DUM")]
    assert len(appearing) == 2


def test_itp_target_missing_file_raises():
    # An ITP target no longer requires a structure (from_itp RDKit-embeds the
    # geometry), but a missing ITP file must still surface clearly.
    from pmx.covalent.mutate import _build_target_residue
    with pytest.raises((FileNotFoundError, OSError)):
        _build_target_residue("SEP", FF, {}, None, None,
                              structure=None, itp="does_not_exist.itp")


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


def _write_tripeptide_internal_ser(tmp_path, FF, strip_h_from=None):
    """Write a three-residue PDB (ALA1-SER2-ALA3) with SER INTERNAL.

    Two constraints force a tripeptide.  A ONE-residue chain is both termini at
    once, which pmx refuses (a dipeptide is the minimum supported input).  And in
    a dipeptide BOTH residues are termini, so a PTM target there would need a
    5-character terminal hybrid name (``S2SPc``), which cannot round-trip through
    a PDB and is refused.  Neither of those is what these tests are about -- they
    exercise the itp-fragment path and the unprotonated-WT warning -- so the
    mutated residue is placed internally.  Geometry only has to be parseable;
    these tests run the in-memory mutate, not pdb2gmx.
    """
    from pmx.covalent import TargetResidue
    out = []
    for resid, rn in ((1, "ALA"), (2, "SER"), (3, "ALA")):
        tr = TargetResidue.from_rtp(rn, FF)
        tmp = str(tmp_path / ("_%s.pdb" % rn))
        tr.write(tmp)
        for l in open(tmp):
            if l[:4] != "ATOM":
                continue
            if strip_h_from == resid and l[76:78].strip() == "H":
                continue
            if strip_h_from == resid and not l[76:78].strip() \
                    and l[12:16].strip().startswith("H"):
                continue
            x = float(l[30:38]) + (resid - 1) * 4.0      # separate the residues
            out.append("%s%4d%s%8.3f%s" % (l[:22], resid, l[26:30], x, l[38:]))
    path = str(tmp_path / "tripeptide.pdb")
    with open(path, "w") as fh:
        fh.writelines(out)
        fh.write("END\n")
    return path


def test_mutate_via_itp_fragment(tmp_path):
    """End-to-end itp path: capped fragment -> hybrid with phosphate dummies."""
    from pmx.gmx import set_gmxlib
    from pmx.model import Model
    from pmx.covalent import TargetResidue, mutate
    set_gmxlib()
    itp, pdb = _write_capped_sep_fragment(tmp_path)
    m = Model(_write_tripeptide_internal_ser(tmp_path, FF))
    name = mutate(m, 2, "SEP", FF, itp=itp, structure=pdb,
                  inplace=True, write_files=False)
    names = [a.name for a in m.residues[1].atoms]
    # SER2 is internal, so this is the ordinary internal hybrid name.
    assert name == "S2SP"
    assert "HN" in names                      # WT backbone H mapped through
    assert any(a.name == "DP" for a in m.residues[1].atoms)  # phosphate appears


def test_protonation_warning_on_hydrogenless_wt(tmp_path):
    from pmx.gmx import set_gmxlib
    from pmx.model import Model
    from pmx.covalent import TargetResidue, mutate
    set_gmxlib()
    # tripeptide with the mutated SER internal, hydrogens stripped from it
    m = Model(_write_tripeptide_internal_ser(tmp_path, FF, strip_h_from=2))
    with pytest.warns(UserWarning, match="unprotonated"):
        mutate(m, 2, "SEP", FF, inplace=True, write_files=False)


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

def test_lone_residue_is_refused(tmp_path):
    """A one-residue chain is both termini; pmx refuses rather than half-patching."""
    from pmx.gmx import set_gmxlib
    from pmx.model import Model
    from pmx.covalent import TargetResidue
    from pmx.alchemy import mutate as alchemy_mutate
    set_gmxlib()
    TargetResidue.from_rtp("SER", FF).write(str(tmp_path / "lone.pdb"))
    m = Model(str(tmp_path / "lone.pdb"))
    with pytest.raises(ValueError, match="dipeptide is the minimum"):
        alchemy_mutate(m, 1, "ALA", FF, inplace=True)
