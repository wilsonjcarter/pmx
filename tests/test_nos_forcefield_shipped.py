"""pmx must ship working NOS / SONOS / sulfenic support, not just the rtp entries.

Two gaps found 2026-09-30, both of which made a sulfenic-cysteine (CSO) system
unbuildable even though every tree already carried the CSO rtp entry:

1. ``charmm36m-mut.ff`` had no CSO entry in its ``.hdb``.  pdb2gmx -ignh then strips
   CSO's hydrogens and cannot rebuild them: *"There were 5 missing atoms in molecule
   Protein"*.  charmm36m-ecc-mut.ff had the entry; only the plain tree was missing it.
2. Neither tree had ``CG321 CT1 C N`` -- the CB-CA-C-N(+1) backbone dihedral for a CSO
   whose next residue is PROLINE (proline's N is type ``N``, not ``NH1``).  Only the
   NH1 case was covered, so a CSO immediately before a Pro failed with *"No default
   Proper Dih. types"*.  Hit on 1puf, CSO18 followed by PRO19.

These are static checks on the shipped data, so they need neither gmx nor a structure.
"""

import glob
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))

MUTFF = os.path.join(HERE, "..", "src", "pmx", "data", "mutff")

_SUB = {"atoms", "bonds", "angles", "dihedrals", "impropers", "cmap",
        "exclusions", "bondedtypes"}


def _rtp_entries(ffdir):
    """-> {entry_name: [atom names]} across every .rtp in the force field."""
    out = {}
    for p in glob.glob(os.path.join(ffdir, "*.rtp")):
        cur = sec = None
        for l in open(p, errors="replace"):
            s = l.split(";")[0].strip()
            if s.startswith("[") and s.endswith("]"):
                nm = s[1:-1].strip()
                if nm in _SUB:
                    sec = nm
                else:
                    cur, sec = nm, None
                continue
            if cur and sec == "atoms" and s:
                f = s.split()
                if len(f) >= 4:
                    out.setdefault(cur, []).append(f[0])
    return out


def _hdb_names(ffdir):
    """-> set of residues with an H-database entry (header is '<res> <count>')."""
    got = set()
    for p in glob.glob(os.path.join(ffdir, "*.hdb")):
        for l in open(p, errors="replace"):
            f = l.split()
            if len(f) == 2 and f[1].isdigit() and not l[0].isspace():
                got.add(f[0])
    return got


def _bonded_keys(ffdir):
    """-> set of (section, (types...)) over ffbonded.itp and everything it includes."""
    arity = {"bondtypes": 2, "angletypes": 3, "dihedraltypes": 4,
             "constrainttypes": 2}
    keys = set()
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
            n = arity.get(sec, 0)
            if not n or not s:
                continue
            f = s.split()
            if len(f) >= n:
                keys.add((sec, tuple(f[:n])))
                keys.add((sec, tuple(reversed(f[:n]))))
    return keys


def _trees_with_cso():
    return [d for d in sorted(glob.glob(os.path.join(MUTFF, "*.ff")))
            if "CSO" in _rtp_entries(d)]


def test_some_tree_ships_cso():
    trees = _trees_with_cso()
    assert trees, "no shipped force field defines CSO -- has the NOS data moved?"


def test_every_tree_with_cso_can_build_its_hydrogens():
    """gap 1: an rtp entry without an .hdb entry cannot be built with -ignh."""
    missing = [os.path.basename(d) for d in _trees_with_cso()
               if "CSO" not in _hdb_names(d)]
    assert not missing, (
        "these force fields define CSO in an .rtp but have no CSO .hdb entry, so "
        "pdb2gmx -ignh cannot rebuild its hydrogens: %s" % missing)


def test_every_tree_with_cso_resolves_its_sulfenic_bonded_terms():
    """The CGenFF-typed terms CSO needs, including the CHARMM/CGenFF junction."""
    need = [("bondtypes", ("CT1", "CG321")),
            ("bondtypes", ("SG301", "OG311")),
            ("angletypes", ("CT1", "CG321", "SG301")),
            ("angletypes", ("CG321", "SG301", "OG311")),
            ("angletypes", ("SG301", "OG311", "HGP1")),
            ("dihedraltypes", ("CT1", "CG321", "SG301", "OG311")),
            ("dihedraltypes", ("HGA2", "CG321", "SG301", "OG311")),
            ("dihedraltypes", ("CG321", "SG301", "OG311", "HGP1"))]
    bad = {}
    for d in _trees_with_cso():
        keys = _bonded_keys(d)
        miss = [" ".join(t) for sec, t in need if (sec, t) not in keys]
        if miss:
            bad[os.path.basename(d)] = miss
    assert not bad, "sulfenic bonded terms unresolved: %s" % bad


def test_cso_before_proline_dihedral_is_present():
    """gap 2: proline's backbone N is type 'N', not 'NH1'."""
    bad = [os.path.basename(d) for d in _trees_with_cso()
           if ("dihedraltypes", ("CG321", "CT1", "C", "N")) not in _bonded_keys(d)]
    assert not bad, (
        "CG321-CT1-C-N missing, so a CSO immediately before a PROLINE cannot be "
        "typed (grompp: 'No default Proper Dih. types'): %s" % bad)


def test_nos_linkage_bonded_terms_are_present():
    """The NOS bond itself: Cys-S-O-N-Lys."""
    need = [("bondtypes", ("NG311", "OG301")),
            ("bondtypes", ("OG301", "SG301")),
            ("angletypes", ("NG311", "OG301", "SG301"))]
    bad = {}
    for d in _trees_with_cso():
        keys = _bonded_keys(d)
        miss = [" ".join(t) for sec, t in need if (sec, t) not in keys]
        if miss:
            bad[os.path.basename(d)] = miss
    assert not bad, "NOS linkage terms unresolved: %s" % bad


def _rtp_bonds(ffdir, name):
    """-> (atom->type map, [(a,b)]) for one rtp entry."""
    import glob
    types = {}; bonds = []
    for p in glob.glob(os.path.join(ffdir, "*.rtp")):
        cur = sec = None
        for l in open(p, errors="replace"):
            s = l.split(";")[0].strip()
            if s.startswith("[") and s.endswith("]"):
                nm = s[1:-1].strip()
                if nm in _SUB:
                    sec = nm
                else:
                    cur, sec = nm, None
                continue
            if cur != name or not s:
                continue
            f = s.split()
            if sec == "atoms" and len(f) >= 4:
                types[f[0]] = f[1]
            elif sec == "bonds" and len(f) >= 2:
                bonds.append((f[0], f[1]))
        if types:
            break
    return types, bonds


CROSSLINK_ENDPOINTS = ("CN0", "KN0", "KS0")


def test_crosslink_endpoints_are_shipped_wherever_nos_is():
    """CN0/KN0/KS0 are the reacted NOS/SONOS endpoints decouple_mutate morphs to.

    They lived only in the PLAIN trees' nos.rtp for a long time, so an ECC cross-link
    build failed with "cross-link target is not in the force field" even though every
    atom type it needed was already defined.  Ported 2026-09-30.
    """
    trees = [d for d in sorted(glob.glob(os.path.join(MUTFF, "*.ff")))
             if os.path.exists(os.path.join(d, "nos.rtp"))]
    assert trees, "no shipped force field carries nos.rtp"
    bad = {}
    for d in trees:
        got = _rtp_entries(d)
        miss = [r for r in CROSSLINK_ENDPOINTS if r not in got]
        if miss:
            bad[os.path.basename(d)] = miss
    assert not bad, "force fields with nos.rtp but missing endpoints: %s" % bad


def test_crosslink_endpoint_atom_types_are_defined():
    for d in sorted(glob.glob(os.path.join(MUTFF, "*.ff"))):
        got = _rtp_entries(d)
        if not any(r in got for r in CROSSLINK_ENDPOINTS):
            continue
        defined = set()
        for p in (glob.glob(os.path.join(d, "*.atp"))
                  + glob.glob(os.path.join(d, "ffnonbonded.itp"))):
            sec = None
            for l in open(p, errors="replace"):
                s = l.split(";")[0].strip()
                if s.startswith("["):
                    sec = s.strip("[] ")
                    continue
                if s and (p.endswith(".atp") or sec == "atomtypes"):
                    defined.add(s.split()[0])
        for r in CROSSLINK_ENDPOINTS:
            if r not in got:
                continue
            types, _ = _rtp_bonds(d, r)
            miss = sorted({t for t in types.values()
                           if t not in defined and not t.startswith("DUM_")})
            assert not miss, "%s/%s references undefined atom types %s" % (
                os.path.basename(d), r, miss)


def test_crosslink_endpoint_internal_bonds_are_typed():
    """Every bond INSIDE CN0/KN0/KS0 must have a bondtype, or grompp cannot build it."""
    bad = {}
    for d in sorted(glob.glob(os.path.join(MUTFF, "*.ff"))):
        got = _rtp_entries(d)
        if not any(r in got for r in CROSSLINK_ENDPOINTS):
            continue
        keys = _bonded_keys(d)
        for r in CROSSLINK_ENDPOINTS:
            if r not in got:
                continue
            types, bonds = _rtp_bonds(d, r)
            miss = []
            for a, b in bonds:
                ta, tb = types.get(a), types.get(b)
                if not ta or not tb:
                    continue          # bond to a neighbouring residue (+N/-C)
                if ("bondtypes", (ta, tb)) not in keys:
                    miss.append("%s-%s (%s-%s)" % (a, b, ta, tb))
            if miss:
                bad["%s/%s" % (os.path.basename(d), r)] = miss
    assert not bad, "untyped internal bonds in the cross-link endpoints: %s" % bad

