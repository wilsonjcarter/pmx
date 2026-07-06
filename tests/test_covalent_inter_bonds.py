#!/usr/bin/env python
"""Tests for the ``[ inter_bonds ]`` MTP section: writer, parser round-trip, and
the ``gen_hybrid_top`` consumer (``_apply_inter_bonds``).

The section serializes only the *forming bond(s)* + mode; ``gen_hybrid_top``
re-resolves the atom names against the parsed topology and derives the induced
angles/dihedrals from the graph.  These tests use hand-built objects (no RDKit /
force field), so they pin the format and the apply logic directly.
"""
import io

from pmx.atom import Atom


# ── writer + parser round-trip ────────────────────────────────────────────────

class _R:
    def __init__(self, atoms, resname, resnA, resnB):
        self.atoms, self.resname, self.resnA, self.resnB = atoms, resname, resnA, resnB


def _atom(name, i, t="CT2", tB="DUM_CT2"):
    a = Atom(name=name, id=i)
    a.atomtype, a.atomtypeB = t, tB
    a.nameB, a.q, a.qB, a.m, a.mB = name, 0.0, 0.0, 12.0, 12.0
    a.x = [0.0, 0.0, 0.0]
    a.type, a.typeB = t, tB
    return a


def test_writer_emits_inter_bonds_section():
    from pmx.nsaa import _make_mtp
    atoms = [_atom("SG2", 1)]
    r = _R(atoms, "C2CY", "CYS", "CYS2")
    txt = "".join(_make_mtp(r, [], [], inter_bonds=[("form", "SG2", "SG2")]))
    assert "[ inter_bonds ]" in txt
    assert "form" in txt and "SG2" in txt and "@SG2" in txt


def test_writer_omits_section_when_no_inter_bonds():
    from pmx.nsaa import _make_mtp
    r = _R([_atom("SG2", 1)], "C2CY", "CYS", "CYS2")
    txt = "".join(_make_mtp(r, [], []))
    assert "[ inter_bonds ]" not in txt


def test_parser_reads_inter_bonds_roundtrip():
    from pmx.nsaa import _make_mtp
    from pmx.mutdb import read_new_mtp_entry
    atoms = [_atom("SG2", 1)]
    r = _R(atoms, "C2CY", "CYS", "CYS2")
    txt = "".join(_make_mtp(r, [], [], inter_bonds=[("form", "SG2", "SG2")]))
    mol, _b, _i, _d, _rot = read_new_mtp_entry("C2CY", filename=io.StringIO(txt))
    assert mol.inter_bonds == [{"mode": "form", "self": "SG2", "partner": "SG2", "cutoff": None}]


# ── consumer: _apply_inter_bonds over a synthetic topology ────────────────────

class _BP:
    def get_bond_param(self, *t):     return (1, 2.05, 300.0)
    def get_angle_param(self, *t):    return (5, 110.0, 50.0, 0.0, 0.0)
    def get_dihedral_param(self, *t): return [(9, 180.0, 2.0, 2)]


class _Topol:
    def __init__(self, bonds):
        self.bonds = [list(b) for b in bonds]
        self.angles, self.dihedrals = [], []
        self.BondedParams = _BP()


class _Res:
    def __init__(self, resname, atoms):
        self.resname, self.atoms = resname, atoms


class _Entry:
    """Stand-in for a parsed MTP entry (only .inter_bonds is read)."""
    def __init__(self, inter_bonds):
        self.inter_bonds = inter_bonds


def _mkres(resname, base_id):
    names = ["CA", "CB2", "SG2", "HB12", "HB22"]
    atoms = {}
    for k, n in enumerate(names):
        a = Atom(name=n, id=base_id + k)
        a.type = a.typeB = "X"
        a.resname = resname
        atoms[n] = a
    return _Res(resname, list(atoms.values())), atoms


def test_apply_inter_bonds_disulfide_form():
    from pmx.alchemy import _apply_inter_bonds

    r1, a1 = _mkres("C2CY", 1)
    r2, a2 = _mkres("C2CY", 100)
    # intra-residue bonds for both sidechains (needed for the graph walk)
    bonds = []
    for a in (a1, a2):
        bonds += [[a["CA"], a["CB2"], 1], [a["CB2"], a["SG2"], 1],
                  [a["CB2"], a["HB12"], 1], [a["CB2"], a["HB22"], 1]]
    topol = _Topol(bonds)

    rlist = [r1, r2]
    rdic = {"C2CY": (_Entry([{"mode": "form", "self": "SG2", "partner": "SG2"}]),
                     [], [], [], {})}

    _apply_inter_bonds(topol, rlist, rdic)

    # exactly one SG2-SG2' bond (symmetric decls deduped)
    assert len(topol.bonds) == len(bonds) + 1
    new_bond = topol.bonds[-1]
    assert {new_bond[0].name, new_bond[1].name} == {"SG2"}
    assert new_bond[3][2] == 0.0 and new_bond[4][2] == 300.0   # A off, B on
    # 2 induced angles, and the graph-complete 7 induced dihedrals
    assert len(topol.angles) == 2
    assert len(topol.dihedrals) == 7


def test_covalent_edit_emit_serializes_inter_bonds_section():
    """The generation side: CovalentEdit.emit() writes an [ inter_bonds ] section
    for the bond it owns, round-tripping through the parser."""
    from pmx.covalent import CovalentEdit
    from pmx.covalent.edit import InterResidueBondMorph, Interaction
    from pmx.mutdb import read_new_mtp_entry

    sg1 = _atom("SG2", 9, "SM", "DUM_SM")
    sg2 = _atom("SG2", 19, "SM", "DUM_SM")
    edit = CovalentEdit()
    # emulate two decouples + a bond on the bound-branch SG2 atoms
    edit.decoupled_render.append(("C2CY", [sg1], "CYS", "CYS2"))
    edit.decoupled_render.append(("C2CY", [sg2], "CYS", "CYS2"))
    edit.bond(InterResidueBondMorph("form", [Interaction("bond", [sg1, sg2])]))

    txt = edit.emit()
    assert txt.count("[ inter_bonds ]") == 1        # declared once, on the owner
    mol, *_ = read_new_mtp_entry("C2CY", filename=io.StringIO(txt))
    assert mol.inter_bonds == [{"mode": "form", "self": "SG2", "partner": "SG2", "cutoff": None}]


def test_apply_inter_bonds_multibond_sonos():
    """Consumer handles a multi-bond cross-link: a Lys declaring two forming
    bonds to two Cys (SONOS-like), including the shared-nitrogen cross-angle."""
    from pmx.alchemy import _apply_inter_bonds

    def lys(base):
        names = ["CE2", "NZ2", "HZ21"]
        d = {}
        for k, n in enumerate(names):
            a = Atom(name=n, id=base + k); a.type = a.typeB = "X"; a.resname = "K2SO"
            d[n] = a
        return _Res("K2SO", list(d.values())), d

    def cys(resname, base, od):
        names = ["CB2", "SG2", od]
        d = {}
        for k, n in enumerate(names):
            a = Atom(name=n, id=base + k); a.type = a.typeB = "X"; a.resname = resname
            d[n] = a
        return _Res(resname, list(d.values())), d

    rl, l = lys(1)
    rca, ca = cys("C2SA", 50, "ODA")
    rcb, cb = cys("C2SB", 80, "ODB")
    bonds = [
        [l["CE2"], l["NZ2"], 1], [l["NZ2"], l["HZ21"], 1],
        [ca["CB2"], ca["SG2"], 1], [ca["SG2"], ca["ODA"], 1],
        [cb["CB2"], cb["SG2"], 1], [cb["SG2"], cb["ODB"], 1],
    ]
    topol = _Topol(bonds)
    rlist = [rl, rca, rcb]
    rdic = {
        "K2SO": (_Entry([{"mode": "form", "self": "NZ2", "partner": "ODA"},
                         {"mode": "form", "self": "NZ2", "partner": "ODB"}]),
                 [], [], [], {}),
        "C2SA": (_Entry([]), [], [], [], {}),
        "C2SB": (_Entry([]), [], [], [], {}),
    }

    _apply_inter_bonds(topol, rlist, rdic)

    assert len(topol.bonds) == len(bonds) + 2          # two forming bonds
    # the ODA-NZ2-ODB cross-angle (both bonds meet at NZ2) is generated
    ang_names = {frozenset((a[0].name, a[2].name)) for a in topol.angles}
    assert frozenset(("ODA", "ODB")) in ang_names


def test_apply_inter_bonds_explicit_gentop_bond():
    """gentop --bond forms a bond deterministically by (resid, atom), even with
    no [ inter_bonds ] section -- the reliable path for far-apart reacting atoms."""
    from pmx.alchemy import _apply_inter_bonds

    r1, a1 = _mkres("C2QN", 1); r1.id = 32
    r2, a2 = _mkres("C2QN", 100); r2.id = 35
    for a in (a1, a2):
        a["CA"].resnr = None
    bonds = []
    for a in (a1, a2):
        bonds += [[a["CA"], a["CB2"], 1], [a["CB2"], a["SG2"], 1]]
    topol = _Topol(bonds)
    rlist = [r1, r2]
    rdic = {"C2QN": (_Entry([]), [], [], [], {})}   # no rule -> only --bond drives it

    n0 = len(topol.bonds)
    # plain atom name "SG" maps to the bound branch "SG2"
    _apply_inter_bonds(topol, rlist, rdic, explicit_bonds=[
        {"resid_a": 32, "atom_a": "SG", "resid_b": 35, "atom_b": "SG",
         "mode": "form"}])
    assert len(topol.bonds) == n0 + 1               # the SG2-SG2 bond was formed
    formed = topol.bonds[-1]
    assert {formed[0].name, formed[1].name} == {"SG2"}


def test_apply_bond_morph_skips_missing_param_with_warning():
    """A generated induced term with no FF parameter is skipped + warned, not a
    crash (surfaces gaps cleanly when params are incomplete/placeholder)."""
    import warnings
    from pmx.covalent.edit import (apply_bond_morph, InterResidueBondMorph,
                                   Interaction)

    class _BPNone:
        def get_bond_param(self, *t):     return (1, 2.0, 300.0)
        def get_angle_param(self, *t):    return None      # missing
        def get_dihedral_param(self, *t): return []

    class _T:
        def __init__(self):
            self.bonds, self.angles, self.dihedrals = [], [], []
            self.BondedParams = _BPNone()

    a = _atom("NZ2", 1, "NG301", "NG301")
    b = _atom("OD2", 2, "OG301", "OG301")
    c = _atom("OD2", 3, "OG301", "OG301")
    morph = InterResidueBondMorph("form", [
        Interaction("bond", [a, b]),
        Interaction("angle", [b, a, c]),          # OD-NZ-OD, no FF param
    ])
    topol = _T()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        apply_bond_morph(topol, morph)
    assert len(topol.bonds) == 1                   # the bond still formed
    assert len(topol.angles) == 0                  # the unparameterised angle skipped
    assert any("SKIPPED" in str(x.message) for x in w)


def test_apply_inter_bonds_noop_without_section():
    from pmx.alchemy import _apply_inter_bonds
    r1, a1 = _mkres("SOME", 1)
    topol = _Topol([[a1["CA"], a1["CB2"], 1]])
    # rdic entry without inter_bonds -> consumer is a no-op
    rdic = {"SOME": (_Entry([]), [], [], [], {})}
    _apply_inter_bonds(topol, [r1], rdic)
    assert topol.angles == [] and topol.dihedrals == [] and len(topol.bonds) == 1
