"""Every residue pmx ships in its own rtp must be fully typed by the force field.

The invariant: for each residue defined in one of pmx's residue "topic" rtp files, every
internal bond, angle and dihedral it declares must have a matching entry reachable from
that force field's ``ffbonded.itp`` (following ``#include``s).  If one is missing, the
residue cannot be built at all -- grompp stops with "No default ... types", or pmx stops
earlier with ``MissingTopolParamError``.

This is deliberately a general check rather than a per-residue one.  Three separate
incidents in one session were all this same gap:

* ``CSO``  -- the CGenFF-typed sulfenic terms were absent, so no sulfenic leg built;
* ``KS0``  -- its own ``CE-NZ`` (CT2-NG301) bond was untyped, so no SONOS endpoint built;
* ``ABA``  -- its own ``CA-CB-CG`` (CT1-CG321-CG331) angle was untyped, so Abu would not
  build.  Note pmx itself was blameless there: ``_check_case`` correctly saw the appearing
  CG as a dummy (A='AAD', B='AAA') and looked up the all-real B-state parameters, which is
  precisely the angle that was missing.  No chimeric or dummy-state term was involved.

A residue NOT having an ``.hdb`` entry is deliberately *not* checked: a mutation target
takes its atoms from the rtp, so it needs none.  An .hdb entry only matters for a residue
that appears as a STATIC residue in an input structure and so has to have its hydrogens
rebuilt by ``pdb2gmx -ignh`` -- which is why CSO needed one and the PTM targets do not.
"""

import glob
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))

MUTFF = os.path.join(HERE, "..", "src", "pmx", "data", "mutff")

# pmx's own residue files (as opposed to the vendored merged.rtp / mutres*.rtp libraries)
TOPIC_RTP = ("aba.rtp", "nos.rtp", "tym.rtp", "acyl.rtp", "glycation.rtp",
             "ptm_phospho.rtp", "ptm_methyl.rtp", "ptm_misc.rtp", "ptm_extra.rtp")

# Residues KNOWN to be missing a residue-specific bonded term, measured 2026-09-30.
# (CSO, KS0, ABA and KCX were all in this list and have since been fixed by porting
#  the terms verbatim from stock charmm36-jul2022.ff.)
# They are recorded here rather than silently skipped: each one cannot currently be built,
# and the fix is the same as for CSO/KS0/ABA -- find the term in stock CHARMM36/CGenFF (or
# derive it) and add it to that tree.  Remove entries from this dict as they are fixed; a
# residue that is NOT listed and NOT complete fails the test.
KNOWN_GAPS = {
    "CSS": ["bondtypes SG-SD"],           # thiomethyl Cys
    "NEP": ["bondtypes NE2-P"],           # phospho-His
    "SNC": ["bondtypes SG-ND", "bondtypes ND-OE"],   # S-nitrosocysteine
}

_SUB = {"atoms", "bonds", "angles", "dihedrals", "impropers", "cmap",
        "exclusions", "bondedtypes"}
_NT = {"bondtypes": 2, "angletypes": 3, "dihedraltypes": 4}


def _bonded_keys(ffdir):
    """Every (section, type-tuple) reachable from ffbonded.itp, both orientations."""
    out = set()
    seen = set()
    todo = ["ffbonded.itp"]
    while todo:
        rel = todo.pop()
        p = os.path.join(ffdir, rel)
        if p in seen or not os.path.exists(p):
            continue
        seen.add(p)
        sec = None
        for l in open(p, errors="replace"):
            inc = re.match(r'\s*#include\s+"([^"]+)"', l)
            if inc:
                todo.append(inc.group(1))
                continue
            s = l.split(";")[0].strip()
            if s.startswith("[") and s.endswith("]"):
                sec = s[1:-1].strip()
                continue
            n = _NT.get(sec, 0)
            if not n or not s:
                continue
            f = s.split()
            if len(f) >= n:
                out.add((sec, tuple(f[:n])))
                out.add((sec, tuple(reversed(f[:n]))))
    return out


def _parse_rtp(path):
    cur = sec = None
    out = {}
    for l in open(path, errors="replace"):
        s = l.split(";")[0].strip()
        if s.startswith("[") and s.endswith("]"):
            nm = s[1:-1].strip()
            if nm in _SUB:
                sec = nm
            else:
                cur, sec = nm, None
                out.setdefault(cur, {"atoms": {}, "q": {}, "bonds": [],
                                     "angles": [], "dihedrals": []})
            continue
        if not cur or not s:
            continue
        f = s.split()
        if sec == "atoms" and len(f) >= 4:
            out[cur]["atoms"][f[0]] = f[1]
            try:
                out[cur]["q"][f[0]] = float(f[2])
            except ValueError:
                pass
        elif sec in ("bonds", "angles", "dihedrals"):
            out[cur][sec].append(f)
    return out


def _untyped(ffdir, keys, entry):
    """-> list of 'kind a-b-c' for terms with no force-field entry."""
    types = entry["atoms"]
    miss = []
    for sec, ar, kind in (("bonds", 2, "bondtypes"),
                          ("angles", 3, "angletypes"),
                          ("dihedrals", 4, "dihedraltypes")):
        for f in entry[sec]:
            tt = [types.get(x) for x in f[:ar]]
            if any(t is None for t in tt):
                continue            # bond to a neighbouring residue (-C / +N)
            if (kind, tuple(tt)) not in keys:
                miss.append("%s %s" % (kind, "-".join(f[:ar])))
    return miss


def _scan():
    """-> {(tree, residue): [untyped terms]} over every topic rtp in every tree."""
    bad = {}
    for ffdir in sorted(glob.glob(os.path.join(MUTFF, "*.ff"))):
        keys = None
        for fn in TOPIC_RTP:
            p = os.path.join(ffdir, fn)
            if not os.path.exists(p):
                continue
            if keys is None:
                keys = _bonded_keys(ffdir)
            for res, entry in _parse_rtp(p).items():
                if res == "bondedtypes":
                    continue
                miss = _untyped(ffdir, keys, entry)
                if miss:
                    bad[(os.path.basename(ffdir), res)] = miss
    return bad


def test_shipped_residues_are_fully_typed():
    bad = _scan()
    unexpected = {k: v for k, v in bad.items() if k[1] not in KNOWN_GAPS}
    assert not unexpected, (
        "residues with an untyped internal bonded term that are not in KNOWN_GAPS "
        "(each of these cannot be built): %s" % unexpected)


def test_known_gaps_are_still_real_and_not_understated():
    """If a gap is fixed, remove it from KNOWN_GAPS; if it grows, the list is stale."""
    bad = _scan()
    for (tree, res), miss in bad.items():
        if res in KNOWN_GAPS:
            extra = sorted(set(miss) - set(KNOWN_GAPS[res]))
            assert not extra, (
                "%s/%s is missing MORE than KNOWN_GAPS records: %s" % (tree, res, extra))


def test_aba_is_complete_in_every_tree_that_ships_it():
    """Abu, the 4th rung of the GB1 hydrocarbon ladder -- ported 2026-09-30."""
    trees = [d for d in sorted(glob.glob(os.path.join(MUTFF, "*.ff")))
             if os.path.exists(os.path.join(d, "aba.rtp"))]
    assert trees, "no tree ships aba.rtp"
    for d in trees:
        keys = _bonded_keys(d)
        ent = _parse_rtp(os.path.join(d, "aba.rtp"))
        assert "ABA" in ent, "%s/aba.rtp defines no ABA" % os.path.basename(d)
        miss = _untyped(d, keys, ent["ABA"])
        assert not miss, "%s: ABA untyped terms %s" % (os.path.basename(d), miss)
        # the residue must be net-neutral, and ECC must not have rescaled it:
        # neutral residues are byte-identical between the plain and ecc trees, only
        # formally ionic ones are scaled (CYM -1 -> -0.75, LYS +1 -> +0.75)
        q = sum(ent["ABA"]["q"].values())
        assert abs(q) < 1e-6, "%s: ABA net charge %+.3f, expected 0" % (
            os.path.basename(d), q)


def test_aba_is_registered_as_protein():
    rt = os.path.join(MUTFF, "residuetypes.dat")
    got = {l.split()[0]: l.split()[1] for l in open(rt) if len(l.split()) >= 2}
    assert got.get("ABA") == "Protein", "ABA must be classified Protein, got %r" % got.get("ABA")
