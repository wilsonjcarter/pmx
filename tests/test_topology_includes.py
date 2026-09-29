"""The water/ion #include statements gentop writes must not depend on the path.

With GMXLIB set, pdb2gmx writes *absolute* include paths into the topology, so
any filter that matches on the whole ``#include`` line also matches the
directory the user happened to run in.  Two such filters existed in
``TopolBase.write_footer``:

    water : ... and "_" not in line     -> a path containing '_'   dropped the
                                           water #include entirely
    ions  : 'ff' in line and 'ion' in line
                                        -> a path containing 'ion' wrote the
                                           water #include a SECOND time, in the
                                           ion slot

Both silently corrupt a solvated topology: grompp then reports
``No such moleculetype SOL``, or SOL as redefined.  Measured on real builds,
a GMXLIB under ``/var/folders/ct/lxwxb4z10_.../T`` lost the include in 17 of 18
golden cases, and a directory named ``fluorination`` duplicated it.
"""

import io
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))

from pmx.forcefield import (TopolBase, _include_path, WATER_ITPS,
                            _is_water_itp, _is_ion_itp)

# tests/conftest.py imports pmx at collection time, before the sys.path line
# above runs, so an installed pmx elsewhere on the machine wins.  Fail loudly
# rather than quietly testing a different checkout.
import pmx.forcefield as _ff
_want = os.path.realpath(os.path.join(HERE, "..", "src", "pmx"))
_got = os.path.realpath(os.path.dirname(_ff.__file__))
if _got != _want:
    raise RuntimeError(
        "pmx was imported from %s, not from this checkout (%s).\n"
        "Run the tests with PYTHONPATH=<repo>/src." % (_got, _want))


FF = "charmm36m-mut.ff"


def _footer(root):
    """The footer pdb2gmx writes for a peptide, rooted at `root`."""
    ff = os.path.join(root, FF)
    return [
        "; Include Position restraint file",
        "#ifdef POSRES",
        '#include "posre.itp"',
        "#endif",
        "",
        "; Include water topology",
        '#include "%s"' % os.path.join(ff, "tip3p.itp"),
        "",
        "#ifdef POSRES_WATER",
        "[ position_restraints ]",
        "1    1       1000       1000       1000",
        "#endif",
        "",
        "; Include topology for ions",
        '#include "%s"' % os.path.join(ff, "ions.itp"),
    ]


def _emit(root, tmp_path):
    top = TopolBase()
    top.footer = _footer(root)
    out = os.path.join(str(tmp_path), "out.top")
    with io.open(out, "w") as fp:
        top.write_footer(fp, what=["water", "posre_water", "ions"])
    text = io.open(out).read()
    return text.count("tip3p.itp"), text.count("ions.itp")


def test_include_path_extracts_the_quoted_file():
    assert _include_path('#include "a/b/tip3p.itp"') == "a/b/tip3p.itp"
    assert _include_path('  #include   "x.itp"  ') == "x.itp"
    assert _include_path("#include <y.itp>") == "y.itp"
    assert _include_path("[ moleculetype ]") == ""
    assert _include_path("; #include \"z.itp\"") == ""


def test_water_and_ion_includes_written_exactly_once(tmp_path):
    """One water #include and one ion #include, whatever the path spells."""
    for root in ("/data/plain",
                 "/data/plain_dir",            # '_' used to drop water
                 "/data/2026_modag/run_1",     # a real project path
                 "/data/fluorination",         # 'ion' used to duplicate water
                 "/var/folders/ct/lxwxb4z10_z8b8m58mq4zrv00000gn/T/tmp1"):
        water, ions = _emit(root, tmp_path)
        assert water == 1, "%s: %d water includes, expected 1" % (root, water)
        assert ions == 1, "%s: %d ion includes, expected 1" % (root, ions)


def test_relative_include_paths_still_work(tmp_path):
    """pdb2gmx without GMXLIB writes '<ff>.ff/tip3p.itp' relative."""
    water, ions = _emit("", tmp_path)
    assert (water, ions) == (1, 1)


def test_amber_ion_file_is_not_mistaken_for_water(tmp_path):
    """Amber names its ion file ``ions_tip3p.itp`` -- it contains 'tip'.

    GROMACS embeds the water model in the ion file name, so any substring test
    for a water model matches the ion file too.  Five of the fourteen shipped
    force fields name it this way; the other nine use a plain ``ions.itp``.
    Matching exact file names separates the two without having to guess.
    """
    ff = "/data/2026_project/amber14sb-q-mut.ff"
    top = TopolBase()
    top.footer = ['#include "%s/tip3p.itp"' % ff,
                  '#include "%s/ions_tip3p.itp"' % ff]
    out = os.path.join(str(tmp_path), "o.top")
    with io.open(out, "w") as fp:
        top.write_footer(fp, what=["water", "ions"])
    lines = [l for l in io.open(out).read().splitlines() if l.strip()]
    assert len(lines) == 2, lines
    assert lines[0].endswith('/tip3p.itp"'), lines[0]
    assert lines[1].endswith('/ions_tip3p.itp"'), lines[1]


def test_non_forcefield_includes_are_not_mistaken_for_water(tmp_path):
    """A user's own itp must not be emitted in the water or ion slot."""
    top = TopolBase()
    top.footer = ['#include "ligand_tip_fragment.itp"',
                  '#include "my_ions_restraints.itp"']
    out = os.path.join(str(tmp_path), "o.top")
    with io.open(out, "w") as fp:
        top.write_footer(fp, what=["water", "ions"])
    assert io.open(out).read().strip() == ""


def test_water_itps_covers_what_pdb2gmx_can_select():
    """WATER_ITPS must cover every model the shipped watermodels.dat files offer.

    That file is what ``pdb2gmx -water`` offers, so it is the authoritative set
    of names that can reach the footer.  Deriving the expectation from the
    force-field data rather than restating the list means a water model added to
    a force field cannot silently fall out of the water slot.

    Coverage, not equality: WATER_ITPS is deliberately a superset.  'tip5pe' is
    selectable in stock GROMACS but in no mutff force field, and gentop also
    reads topologies built against stock force fields.
    """
    mutff = os.path.join(HERE, "..", "src", "pmx", "data", "mutff")
    selectable = set()
    found = 0
    for ff in sorted(os.listdir(mutff)):
        dat = os.path.join(mutff, ff, "watermodels.dat")
        if not os.path.isfile(dat):
            continue
        found += 1
        for line in io.open(dat):
            line = line.strip()
            if line and not line.startswith(";"):
                selectable.add(line.split()[0])
    assert found > 0, "no watermodels.dat found under %s" % mutff
    uncovered = {n + ".itp" for n in selectable} - set(WATER_ITPS)
    assert not uncovered, (
        "these force fields offer water models WATER_ITPS does not cover: %s"
        % sorted(uncovered))


def test_no_name_is_both_water_and_ion():
    """'ions_tip3p.itp' must land in exactly one slot."""
    for model in sorted(n[:-4] for n in WATER_ITPS):
        for name in (model + ".itp", "ions_" + model + ".itp"):
            assert not (_is_water_itp(name) and _is_ion_itp(name)), name
        assert _is_water_itp(model + ".itp")
        assert _is_ion_itp("ions_" + model + ".itp")
    assert _is_ion_itp("ions.itp")
    assert not _is_water_itp("ions.itp")


def test_an_include_is_added_only_if_missing(tmp_path):
    """The same file included twice in the footer goes out once."""
    ff = "/data/x.ff"
    top = TopolBase()
    top.footer = ['#include "%s/tip3p.itp"' % ff,
                  '#include "%s/tip3p.itp"' % ff,
                  '#include "%s/ions_tip3p.itp"' % ff]
    out = os.path.join(str(tmp_path), "o.top")
    with io.open(out, "w") as fp:
        top.write_footer(fp, what=["water", "ions"])
    text = io.open(out).read()
    assert text.count("/tip3p.itp") == 1, text
    assert text.count("ions_tip3p.itp") == 1, text


def test_write_footer_all_is_not_broken(tmp_path):
    """what=['all'] is the default signature; it used to raise NameError."""
    top = TopolBase()
    top.footer = ['; c', '#ifdef POSRES', '#include "posre.itp"', '#endif',
                  '#include "x.ff/tip3p.itp"', '#endif']
    out = os.path.join(str(tmp_path), "o.top")
    with io.open(out, "w") as fp:
        top.write_footer(fp, what=["all"])
    text = io.open(out).read()
    assert "tip3p.itp" in text
    # repeated non-include lines are structure, not duplicates -- keep them
    assert text.count("#endif") == 2, text
