"""A reverse hybrid must be indistinguishable from doing that mutation forward.

Reversing X->Y *is* the mutation Y->X, so it must come out with the same name and
the same atoms as asking for Y->X directly.  Two bugs followed from it not doing so:

1. The name front-loaded the state-A residue and truncated to four characters,
   ``(target + '2' + source_letter)[:4]``.  Every canonical residue name is three
   letters, so the truncation always chopped off exactly the character identifying
   the source: EVERY reverse morph into ALA was ``ALA2``, so ALA->CYS and ALA->LYS
   collided, and the name carried no information about where it came from.
2. Disappearing hydrogens were renamed ``('D' + name)[:4]``, which is not injective
   for a 4-character name: ARG's HH11/HH12 both became DHH1.  pdb2gmx kept one atom
   of each pair and state B silently lost 0.92 e.

Both vanish once the reverse path spells things the way the shipped library does:
hydrogen dummies take pmx's ``HV<n>`` counter (what the forward generator and the
pre-built hybrids both use), and the name becomes the library's name for the inverse
morph.  All 380 canonical ``X2Y`` entries in mutres.rtp have their inverse ``Y2X``
shipped, so the reverse hybrid lines up with the curated entry instead of shadowing
it under a second name.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))

import pmx.library as _L

_want = os.path.realpath(os.path.join(HERE, "..", "src", "pmx"))
_got = os.path.realpath(os.path.dirname(_L.__file__))
if _got != _want:
    raise RuntimeError("tests imported pmx from %s, expected %s" % (_got, _want))

from pmx.library import hybrid_name, _canonical_code, _one_letter


# reverse morph -> the name the shipped library uses for the inverse mutation
CANONICAL = {
    ("LYS", "LEU"): "L2K",
    ("ARG", "ALA"): "A2R",
    ("THR", "VAL"): "V2T",
    ("CYS", "ALA"): "A2C",
    ("LYS", "ALA"): "A2K",
    ("MET", "LYS"): "K2M",
    ("TRP", "GLY"): "G2W",
}


def test_reverse_of_a_canonical_morph_is_the_library_name():
    for (src, tgt), want in sorted(CANONICAL.items()):
        got = hybrid_name(src, tgt, reverse=True)
        assert got == want, "%s->%s reverse: got %r, want %r" % (src, tgt, got, want)


def test_reverse_is_not_the_old_front_loaded_form():
    """The regression that produced ALA2 / LEU2 / VAL2."""
    for (src, tgt) in CANONICAL:
        assert hybrid_name(src, tgt, reverse=True) != tgt + "2"


def test_two_reverse_morphs_into_the_same_residue_no_longer_collide():
    """C41A and K59A both needed a reverse into ALA; both used to be 'ALA2'."""
    a = hybrid_name("CYS", "ALA", reverse=True)
    b = hybrid_name("LYS", "ALA", reverse=True)
    assert a != b, (a, b)


def test_no_collisions_across_the_whole_canonical_set():
    """No reverse name may denote two different morphs.

    Keyed on the canonical CODE pair, not the residue names: CYS/CYN/CYX all share
    the code 'C', so they are the same morph in library terms and SHOULD share a name.
    """
    seen = {}
    res = [r for r in _one_letter if len(r) == 3]
    for src in res:
        for tgt in res:
            if src == tgt:
                continue
            cs, ct = _canonical_code(src), _canonical_code(tgt)
            if not cs or not ct:
                continue
            name = hybrid_name(src, tgt, reverse=True)
            morph = (ct, cs)
            if name in seen:
                assert seen[name] == morph, (
                    "reverse name %r denotes both %s and %s" % (name, seen[name], morph))
            else:
                seen[name] = morph
    assert len(seen) > 300, len(seen)


def test_reverse_names_fit_gromacs_four_characters():
    res = [r for r in _one_letter if len(r) == 3]
    for src in res:
        for tgt in res:
            if src == tgt:
                continue
            assert len(hybrid_name(src, tgt, reverse=True)) <= 4


def test_phosphate_protonation_naming_is_untouched():
    """The YY12/YY21 spelling must survive: it is checked by its own test module too,
    but the reverse rule runs near it and must not capture these."""
    assert hybrid_name("YP1", "YP2") == "YY12"
    assert hybrid_name("YP1", "YP2", reverse=True) == "YY21"
    assert hybrid_name("YP2", "YP1", reverse=True) == "YY12"


def test_ptm_target_keeps_the_front_loaded_form():
    """A PTM/NSAA endpoint has no classic code, so the inverse has no shipped hybrid
    and delegating would collide: MLZ->LYS would be 'M2LS', which is already MET->LYS."""
    assert hybrid_name("LYS", "MLZ", reverse=True) == "MLZ2"
    assert hybrid_name("MET", "LYS") != hybrid_name("LYS", "MLZ", reverse=True)
