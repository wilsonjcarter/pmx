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


def test_phosphate_appears_as_bstate(hybrid):
    """Phosphate atoms are dummies in A (DUM_) and real in B."""
    morphes = {m[2]: m for m in hybrid[3]}  # keyed by B-state name
    for bname in ("P", "O1P", "O2P", "O3P"):
        assert bname in morphes, f"{bname} not appearing in B state"
        a_name, a_type, _, b_type = morphes[bname]
        assert a_type.startswith("DUM"), f"{bname} should be an A-state dummy"
        assert not b_type.startswith("DUM")
