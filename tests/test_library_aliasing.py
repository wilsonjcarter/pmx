#!/usr/bin/env python
import copy

import numpy as np

from pmx import builder, library
from pmx.chain import Chain
from pmx.molecule import Molecule


def _library_unchanged(snap):
    for k in snap:
        for e, s in zip(library._aacids[k], snap[k]):
            if not np.allclose(e[1], s[1]):
                return False
    return True


def test_new_aa_does_not_modify_library():
    snap = copy.deepcopy(library._aacids)
    r = Molecule().new_aa('TRP')
    r.translate([10., 0., 0.])
    assert _library_unchanged(snap)


def test_make_residue_does_not_modify_library():
    snap = copy.deepcopy(library._aacids)
    r = builder.make_residue('TRP')
    for a in r.atoms:
        a.x[0] += 10.
    assert _library_unchanged(snap)


def test_chain_create_repeatable_after_com():
    def capped(seq):
        c = Chain().create(seq)
        c.add_nterm_cap()
        c.add_cterm_cap()
        return c
    ref = np.array([a.x for a in capped('GGG').atoms], dtype=float)
    for seq in ('GPG', 'GQG', 'GWG'):
        capped(seq).com()
    again = np.array([a.x for a in capped('GGG').atoms], dtype=float)
    np.testing.assert_allclose(again, ref)
