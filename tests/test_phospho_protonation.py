"""Phosphate protonation changes must get their own hybrid names.

A phosphoamino acid has pKa2 ~5.7-6.2, so both the monoanion (-1) and the
dianion (-2) are populated near neutral pH and morphing between them is a
routine calculation.  mutff ships both states for all three phosphorylated
residues (SP1/SP2, THP1/THP2, YP1/YP2).

The generic naming rule cannot express such a morph: ``_target_code('YP2')`` is
``'Y2'``, so ``YP1->YP2`` would be named ``Y2Y2`` -- which is already
``TYR->YP2``.  pmx's cross-morph guard then refuses the build (correctly), which
is how this was found.  ``<X><X><from><to>`` is used instead.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))

from pmx import library as L

import pmx.library as _L
_want = os.path.realpath(os.path.join(HERE, "..", "src", "pmx"))
_got = os.path.realpath(os.path.dirname(_L.__file__))
if _got != _want:
    raise RuntimeError("pmx imported from %s, not this checkout (%s); "
                       "run with PYTHONPATH=<repo>/src" % (_got, _want))

PAIRS = [("SP1", "SP2", "SS"), ("THP1", "THP2", "TT"), ("YP1", "YP2", "YY")]


def test_protonation_pairs_are_detected():
    for a, b, _ in PAIRS:
        assert L.phospho_protonation_pair(a, b) is not None, (a, b)
        assert L.phospho_protonation_pair(b, a) is not None, (b, a)


def test_non_pairs_fall_through():
    # different parent residue, and same state, are both not protonation changes
    assert L.phospho_protonation_pair("YP1", "SP2") is None
    assert L.phospho_protonation_pair("YP1", "YP1") is None
    assert L.phospho_protonation_pair("TYR", "YP2") is None
    assert L.phospho_protonation_pair("LEU", "PHE") is None


def test_names_are_direction_explicit():
    for a, b, stem in PAIRS:
        assert L.hybrid_name(a, b) == stem + "12"
        assert L.hybrid_name(a, b, reverse=True) == stem + "21"
        # and the opposite morph is the mirror image
        assert L.hybrid_name(b, a) == stem + "21"
        assert L.hybrid_name(b, a, reverse=True) == stem + "12"


def test_name_depends_only_on_the_state_pair():
    """The same physical hybrid must always get the same name.

    ``hybrid_name(YP1, YP2, reverse=True)`` and ``hybrid_name(YP2, YP1)`` describe
    the SAME hybrid -- state A is YP2 and state B is YP1 either way -- so they
    must agree.  That canonicality is what stops the accumulator being asked to
    store one morph under two names.
    """
    for a, b, _ in PAIRS:
        assert L.hybrid_name(a, b, reverse=True) == L.hybrid_name(b, a)
        assert L.hybrid_name(b, a, reverse=True) == L.hybrid_name(a, b)


def test_the_two_directions_are_distinct():
    """A -> B and B -> A are different hybrids and must NOT share a name."""
    for a, b, _ in PAIRS:
        assert L.hybrid_name(a, b) != L.hybrid_name(b, a)


def test_does_not_shadow_the_generic_rule():
    """The names that made this necessary must still mean what they meant."""
    assert L.hybrid_name("TYR", "YP1") == "Y2Y1"
    assert L.hybrid_name("TYR", "YP2") == "Y2Y2"
    # The reverse convention changed 2026-09-30: a reverse hybrid is now named the
    # way the shipped library names that mutation, using the canonical codes.  The
    # library spells TYR->YP1 as 'Y2P1' (Y = TYR, P1 = YP1), so the inverse YP1->TYR
    # is 'P12Y'.  The old value here, 'YP12', was the front-loaded form
    # (state-A residue + '2'), which dropped the target and collided -- every reverse
    # morph into ALA was 'ALA2'.  See tests/test_reverse_hybrid_equivalence.py.
    assert L.hybrid_name("TYR", "YP1", reverse=True) == "P12Y"
    assert L.hybrid_name("LEU", "PHE") == "L2PE"
    assert L.hybrid_name("LYS", "MLZ") == "K2MZ"
    # a protonation name must not be produced by any of those
    generic = {L.hybrid_name(s, t, reverse=r)
               for s, t in (("TYR", "YP1"), ("TYR", "YP2"), ("LEU", "PHE"))
               for r in (False, True)}
    assert generic.isdisjoint({"YY12", "YY21", "SS12", "SS21", "TT12", "TT21"})


def test_all_three_residues_have_both_states_in_the_force_field():
    """The naming is only useful if mutff actually ships both states."""
    import glob
    import re
    ff = os.path.join(HERE, "..", "src", "pmx", "data", "mutff",
                      "charmm36m-mut.ff")
    names = set()
    for rtp in glob.glob(os.path.join(ff, "*.rtp")):
        for line in open(rtp, errors="replace"):
            m = re.match(r"^\[\s*([A-Z0-9]{2,4})\s*\]", line)
            if m:
                names.add(m.group(1))
    for a, b, _ in PAIRS:
        assert a in names, "%s missing from %s" % (a, ff)
        assert b in names, "%s missing from %s" % (b, ff)
