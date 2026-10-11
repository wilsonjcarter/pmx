#!/usr/bin/env python
"""Tests for the PTM shorthand resolver (pmx.covalent.ptm)."""
import pytest

from pmx.covalent import ptm


def test_phospho_protonation_states():
    assert ptm.resolve("SER", "p1") == ("SP1", -1)
    assert ptm.resolve("SER", "p2") == ("SP2", -2)
    # THP1/THP2, NOT TP1/TP2: CHARMM's TP1/TP2 are phospho-TYROSINE under
    # different atom naming, so 'T p1' used to install an aromatic ring.
    assert ptm.resolve("THR", "p1") == ("THP1", -1)
    assert ptm.resolve("THR", "p2") == ("THP2", -2)
    assert ptm.resolve("TYR", "p2") == ("YP2", -2)


@pytest.mark.parametrize("wt,code", [("SER", "p1"), ("SER", "p2"),
                                     ("THR", "p1"), ("THR", "p2"),
                                     ("TYR", "p1"), ("TYR", "p2")])
def test_phospho_target_is_the_right_chemistry(wt, code):
    """Phosphorylation adds a phosphate (P + 3 O) and NOTHING else.

    Asserting the table against itself cannot catch a mis-wired target -- that is
    how 'THR p1' came to point at TP1, a phospho-tyrosine whose ring adds eleven
    heavy atoms.  Compare the built residues instead.
    """
    pytest.importorskip("rdkit")
    from pmx.gmx import set_gmxlib
    from pmx.covalent import TargetResidue
    set_gmxlib()
    target, _ = ptm.resolve(wt, code)
    heavy = lambda rn: sum(
        1 for a in TargetResidue.from_rtp(rn, "charmm36m-mut").residue.atoms
        if a.m > 2.0)
    assert heavy(target) - heavy(wt) == 4, (
        f"{wt} + {code} -> {target} changes the heavy-atom count by "
        f"{heavy(target) - heavy(wt)}, not 4 (P + 3 O); the target is not a "
        f"phosphorylated {wt}")


def test_methylation_degrees():
    assert ptm.resolve("LYS", "me1") == ("MLZ", 0)
    assert ptm.resolve("LYS", "me2") == ("MLY", 0)
    assert ptm.resolve("LYS", "me3") == ("M3L", 0)
    # the charged lysine is renamed LYP (amber building block) before the PTM
    # code is resolved -- it is the methylation substrate, so it must resolve too
    assert ptm.resolve("LYP", "me1") == ("MLZ", 0)


def test_case_insensitive():
    assert ptm.resolve("ser", "P2") == ("SP2", -2)
    assert ptm.resolve("Lys", "ME3") == ("M3L", 0)


def test_digit_required_no_proline_collision():
    # bare 'p' / 'me' are NOT PTM codes -> treated as literal residue names
    assert not ptm.is_ptm_code("p")
    assert not ptm.is_ptm_code("me")
    assert not ptm.is_ptm_code("P")        # Proline one-letter code, untouched
    assert ptm.resolve("SER", "p") is None
    assert ptm.resolve("ALA", "PRO") is None
    assert ptm.is_ptm_code("p1") and ptm.is_ptm_code("me2")


def test_non_ptm_token_returns_none():
    assert ptm.resolve("SER", "SEP") is None    # literal target name passes through
    assert ptm.resolve("LYS", "ALA") is None


def test_wrong_residue_for_code_raises():
    with pytest.raises(ValueError, match="phosphorylation"):
        ptm.resolve("LYS", "p2")            # can't phosphorylate a lysine
    with pytest.raises(ValueError, match="methylation"):
        ptm.resolve("SER", "me1")           # no methyl-serine here
    assert ptm.resolve("SER", "p4") is None  # p4 isn't a code -> passthrough
