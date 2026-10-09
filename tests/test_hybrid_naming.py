#!/usr/bin/env python
"""Hybrid residue naming: uniqueness + the anti-collision guard.

Regression for the K2ML defect -- LYS->MLZ and LYS->MLY (and every phospho
mono/di pair) truncated to the same 4-char name and silently overwrote each
other in the shared nsaa.rtp.  These are pure-logic tests (no RDKit / GROMACS).
"""
import os

import importlib

from pmx import library

# NB: import the submodule, not the re-exported ``mutate`` function of the same
# name in pmx.covalent.__init__.
_mutate_mod = importlib.import_module("pmx.covalent.mutate")

# The shipped PTM / covalent target set (source resname, target resname).
PTM_TARGETS = [
    ("LYS", "MLZ"), ("LYS", "MLY"), ("LYS", "M3L"),
    ("SER", "SP1"), ("SER", "SP2"), ("SER", "SEP"),
    ("THR", "TP1"), ("THR", "TP2"),
    ("TYR", "YP1"), ("TYR", "YP2"), ("TYR", "TYM"),
    ("ALA", "GACE"),
]


def test_names_stay_within_the_four_char_cap():
    for src, tgt in PTM_TARGETS:
        for rev in (False, True):
            assert len(library.hybrid_name(src, tgt, reverse=rev)) <= 4


def test_forward_names_unique_across_ptm_set():
    names = [library.hybrid_name(s, t) for s, t in PTM_TARGETS]
    assert len(set(names)) == len(names), dict(zip(
        [f"{s}->{t}" for s, t in PTM_TARGETS], names))


def test_forward_and_reverse_names_never_collide():
    names = ([library.hybrid_name(s, t) for s, t in PTM_TARGETS]
             + [library.hybrid_name(s, t, reverse=True) for s, t in PTM_TARGETS])
    assert len(set(names)) == len(names)


def test_discriminating_char_preserved():
    # the exact pairs the old truncation collapsed
    assert library.hybrid_name("LYS", "MLZ") != library.hybrid_name("LYS", "MLY")
    assert library.hybrid_name("SER", "SP1") != library.hybrid_name("SER", "SP2")
    assert library.hybrid_name("TYR", "YP1") != library.hybrid_name("TYR", "YP2")


def test_accumulate_allows_rebuild_of_same_morph(tmp_path):
    ff = str(tmp_path)
    entry = "[ K2MZ ] ; LYS -> MLZ\n[ atoms ]\n"
    _mutate_mod._accumulate_hybrid(ff, "K2MZ", entry, entry, False)
    # rebuilding the identical morph replaces its own entry -- no error
    _mutate_mod._accumulate_hybrid(ff, "K2MZ", entry, entry, False)
    assert os.path.isfile(os.path.join(ff, "nsaa.rtp"))


def test_accumulate_rejects_cross_morph_collision(tmp_path):
    ff = str(tmp_path)
    _mutate_mod._accumulate_hybrid(
        ff, "K2MZ", "[ K2MZ ] ; LYS -> MLZ\n[ atoms ]\n",
        "[ K2MZ ] ; LYS -> MLZ\n", False)
    # a DIFFERENT morph hashing to the same name must be refused, not overwrite
    try:
        _mutate_mod._accumulate_hybrid(
            ff, "K2MZ", "[ K2MZ ] ; LYS -> MLY\n[ atoms ]\n",
            "[ K2MZ ] ; LYS -> MLY\n", False)
    except ValueError as e:
        assert "K2MZ" in str(e)
    else:
        raise AssertionError("cross-morph collision was not caught")
