"""Reverse-direction dummy names must stay unique within the residue.

Generating a *reverse* hybrid (product as state A, reactant as state B) is
pmx-fork-only: upstream ships hand-curated forward hybrids and never has to
invent a dummy name for a disappearing atom.  The reverse path therefore has to
name those atoms itself, and GROMACS caps an atom name at 4 characters, so
``'D' + name`` has to lose a character for a 4-character name.

Plain truncation is not injective.  ARG's ``HH11``/``HH12`` both truncate to
``DHH1`` (and LEU's ``HD11``/``HD12``/``HD13`` to ``DHD1``), pdb2gmx then keeps
one atom of the pair, and state B silently loses atoms *and their charge* -- an
ALA->ARG reverse hybrid came out at +0.08 e instead of +1.00 e, which also
breaks charge conservation of a double-system single-box leg.

The forward generator already guards this (``nsaa._merge_molecules`` numbers its
hydrogen dummies ``HV1``, ``HV2``, ...); this pins the same property on the
reverse path.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))

import pmx.library as _L

_want = os.path.realpath(os.path.join(HERE, "..", "src", "pmx"))
_got = os.path.realpath(os.path.dirname(_L.__file__))
if _got != _want:
    raise RuntimeError(
        "tests imported pmx from %s, expected %s -- set "
        "PYTHONPATH=%s" % (_got, _want, os.path.join(HERE, "..", "src")))

from pmx.covalent.mutate import _unique_dummy_name


# sidechain atom names of the residues whose 4-character names share a
# 3-character prefix -- the ones plain truncation collapses.
COLLIDING = {
    "ARG": ["CG", "HG1", "HG2", "CD", "HD1", "HD2", "NE", "HE", "CZ",
            "NH1", "HH11", "HH12", "NH2", "HH21", "HH22"],
    "LEU": ["CG", "HG", "CD1", "HD11", "HD12", "HD13",
            "CD2", "HD21", "HD22", "HD23"],
    "VAL": ["CG1", "HG11", "HG12", "HG13", "CG2", "HG21", "HG22", "HG23"],
    "ILE": ["CG1", "HG11", "HG12", "CG2", "HG21", "HG22", "HG23",
            "CD", "HD1", "HD2", "HD3"],
    "THR": ["OG1", "HG1", "CG2", "HG21", "HG22", "HG23"],
    "ASN": ["CG", "OD1", "ND2", "HD21", "HD22"],
    "GLN": ["CG", "HG1", "HG2", "CD", "OE1", "NE2", "HE21", "HE22"],
}


def _assign(name, used):
    out = _unique_dummy_name(name, used)
    used.add(out)
    return out


def test_no_collisions_for_any_residue():
    for resname, names in sorted(COLLIDING.items()):
        used = set()
        out = [_assign(n, used) for n in names]
        assert len(set(out)) == len(names), \
            "%s: %d names collapsed to %d: %s" % (
                resname, len(names), len(set(out)), list(zip(names, out)))


def test_names_fit_gromacs_four_characters():
    for resname, names in sorted(COLLIDING.items()):
        used = set()
        for n in names:
            out = _assign(n, used)
            assert 1 <= len(out) <= 4, "%s: %r -> %r is not 1-4 chars" % (
                resname, n, out)


def test_arg_guanidinium_hydrogens_stay_four_distinct_atoms():
    """The exact case that lost 0.92 e from an ALA->ARG reverse hybrid."""
    used = set()
    out = [_assign(n, used) for n in ("HH11", "HH12", "HH21", "HH22")]
    assert len(set(out)) == 4, out


def test_first_choice_is_unchanged_when_free():
    """The fix must be additive: a name that did not collide keeps its spelling.

    Existing correct hybrids in the shipped force fields were generated with the
    old rule, so anything that was already unique has to come out identical or
    the golden fixtures drift.
    """
    for name in ("CG", "HG1", "NE", "HE", "CZ", "NH1", "NZ", "OG1", "SG"):
        assert _unique_dummy_name(name, set()) == ("D" + name)[:4]


def test_avoids_a_real_atom_name_already_in_the_residue():
    assert _unique_dummy_name("HH11", {"DHH1"}) != "DHH1"


def test_raises_rather_than_returning_a_duplicate():
    """Better a loud failure than a silently dropped atom."""
    import pytest
    # exhaust every spelling the namer can produce for this atom
    used = set()
    while True:
        try:
            used.add(_unique_dummy_name("HH11", used))
        except ValueError:
            break
    # 'DHH1' + 'DHH2'..'DHH9' + 'DH10'..'DH99'  (n=1 re-proposes 'DHH1')
    assert len(used) == 99, sorted(used)
    with pytest.raises(ValueError):
        _unique_dummy_name("HH11", used)
