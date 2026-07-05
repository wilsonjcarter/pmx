#!/usr/bin/env python
"""Tests for the generic, graph-derived inter-residue bond-morph.

``induced_bond_morph`` is the unified replacement for the hand-authored
``disulfide_morph`` / ``nos_morph``: a covalent modification is specified by
*which bonds form*, and the induced angles/dihedrals follow from the molecular
graph -- exactly as GROMACS generates bonded terms from the bond list.

These tests use hand-built graphs (no RDKit / force field needed) so they pin
the enumeration rule directly.
"""
from pmx.atom import Atom
from pmx.covalent.edit import induced_bond_morph, disulfide_morph


def _mk(name, i):
    a = Atom(name=name, id=i)
    a.type = a.typeB = "X"
    return a


def _graph(edges):
    atoms, adj = {}, {}
    for a, b in edges:
        for n in (a, b):
            if n not in atoms:
                atoms[n] = _mk(n, len(atoms) + 1)
    for a, b in edges:
        adj.setdefault(atoms[a], set()).add(atoms[b])
        adj.setdefault(atoms[b], set()).add(atoms[a])
    return atoms, adj


def _names(it):
    return tuple(a.name for a in it.atoms)


def _canon_ang(t):
    return tuple(sorted((t[0], t[2]))) + (t[1],)


def _canon_dih(t):
    return min(t, tuple(reversed(t)))


# bound-branch disulfide graph: two Cys sidechains joined by SG2-SG2'
_DISULFIDE_EDGES = [
    ("CA", "CB2"), ("CB2", "SG2"), ("CB2", "HB12"), ("CB2", "HB22"),
    ("CA'", "CB2'"), ("CB2'", "SG2'"), ("CB2'", "HB12'"), ("CB2'", "HB22'"),
    ("SG2", "SG2'"),
]


def test_induced_bond_and_angles_match_handlist():
    atoms, adj = _graph(_DISULFIDE_EDGES)
    m = induced_bond_morph([(atoms["SG2"], atoms["SG2'"])], adj, mode="form")

    bonds = [_names(i) for i in m.interactions if i.kind == "bond"]
    angs = {_canon_ang(_names(i)) for i in m.interactions if i.kind == "angle"}

    assert bonds == [("SG2", "SG2'")]

    c1 = {"CA": atoms["CA"], "CB2": atoms["CB2"], "SG2": atoms["SG2"]}
    c2 = {"CA": atoms["CA'"], "CB2": atoms["CB2'"], "SG2": atoms["SG2'"]}
    hand = disulfide_morph(c1, c2, mode="form", ca="CA", cb="CB2", sg="SG2")
    hand_ang = {_canon_ang(_names(i)) for i in hand.interactions if i.kind == "angle"}
    assert angs == hand_ang


def test_induced_dihedrals_are_handlist_plus_h_torsions():
    """Graph derivation = hand-list dihedrals + the 4 HB-CB-SG-SG' torsions the
    hand-list omits (which a real GROMACS disulfide does include)."""
    atoms, adj = _graph(_DISULFIDE_EDGES)
    m = induced_bond_morph([(atoms["SG2"], atoms["SG2'"])], adj, mode="form")
    dihs = {_canon_dih(_names(i)) for i in m.interactions if i.kind == "dihedral"}

    c1 = {"CA": atoms["CA"], "CB2": atoms["CB2"], "SG2": atoms["SG2"]}
    c2 = {"CA": atoms["CA'"], "CB2": atoms["CB2'"], "SG2": atoms["SG2'"]}
    hand = disulfide_morph(c1, c2, mode="form", ca="CA", cb="CB2", sg="SG2")
    hand_dih = {_canon_dih(_names(i)) for i in hand.interactions if i.kind == "dihedral"}

    assert hand_dih <= dihs                       # every hand dihedral present
    extra = dihs - hand_dih
    assert extra == {
        ("HB12", "CB2", "SG2", "SG2'"),
        ("HB22", "CB2", "SG2", "SG2'"),
        ("HB12'", "CB2'", "SG2'", "SG2"),
        ("HB22'", "CB2'", "SG2'", "SG2"),
    }


def test_induced_is_deterministic():
    atoms, adj = _graph(_DISULFIDE_EDGES)
    a = induced_bond_morph([(atoms["SG2"], atoms["SG2'"])], adj, mode="form")
    b = induced_bond_morph([(atoms["SG2"], atoms["SG2'"])], adj, mode="form")
    assert [(_names(i), i.kind) for i in a.interactions] == \
           [(_names(i), i.kind) for i in b.interactions]


def test_nos_single_bond_reproduces_nos_morph():
    """A NOS bridge is one inter-residue bond (NZ2-OD2); the whole nos_morph
    interaction set falls out of the graph from that single edge."""
    from pmx.covalent.edit import nos_morph
    edges = [
        # Lys bound branch
        ("CD", "CE2"), ("CE2", "NZ2"), ("CE2", "HE21"), ("CE2", "HE22"),
        ("NZ2", "HZ21"),
        # Cys bound branch with bridging OD2
        ("CA", "CB2"), ("CB2", "SG2"), ("SG2", "OD2"),
        # the forming inter-residue bond
        ("NZ2", "OD2"),
    ]
    atoms, adj = _graph(edges)
    m = induced_bond_morph([(atoms["NZ2"], atoms["OD2"])], adj, mode="form")
    got_ang = {_canon_ang(_names(i)) for i in m.interactions if i.kind == "angle"}
    got_dih = {_canon_dih(_names(i)) for i in m.interactions if i.kind == "dihedral"}

    cys = {"CA": atoms["CA"], "CB2": atoms["CB2"], "SG2": atoms["SG2"],
           "OD2": atoms["OD2"]}
    lys = {"CD": atoms["CD"], "CE2": atoms["CE2"], "NZ2": atoms["NZ2"],
           "HZ21": atoms["HZ21"], "HE21": atoms["HE21"], "HE22": atoms["HE22"]}
    hand = nos_morph(cys, lys, mode="form")
    hand_ang = {_canon_ang(_names(i)) for i in hand.interactions if i.kind == "angle"}
    hand_dih = {_canon_dih(_names(i)) for i in hand.interactions if i.kind == "dihedral"}

    # NOS: the graph reproduces the hand-list exactly (all heavy+H already listed)
    assert got_ang == hand_ang
    assert got_dih == hand_dih


def test_sonos_two_bonds_generate_shared_atom_cross_angle():
    """SONOS bridge (Cys1-S-O-N-O-S-Cys2): two forming bonds meet at the Lys N.
    The graph rule generates the ODa-NZ-ODb cross-angle -- the term a per-bond
    hand-list would miss -- from the two bonds sharing NZ."""
    edges = [
        ("CE", "NZ"), ("NZ", "HZ"),
        ("CBa", "SGa"), ("SGa", "ODa"),
        ("CBb", "SGb"), ("SGb", "ODb"),
        ("NZ", "ODa"), ("NZ", "ODb"),         # the two forming bonds
    ]
    atoms, adj = _graph(edges)
    m = induced_bond_morph(
        [(atoms["NZ"], atoms["ODa"]), (atoms["NZ"], atoms["ODb"])], adj, "form")
    bonds = {_names(i) for i in m.interactions if i.kind == "bond"}
    angs = {_canon_ang(_names(i)) for i in m.interactions if i.kind == "angle"}

    assert bonds == {("NZ", "ODa"), ("NZ", "ODb")}
    # the cross-angle where both forming bonds meet
    assert _canon_ang(("ODa", "NZ", "ODb")) in angs
    # and the per-bond incident angles
    assert _canon_ang(("NZ", "ODa", "SGa")) in angs
    assert _canon_ang(("CE", "NZ", "ODa")) in angs


def test_btk_acrylamide_cysteine_boundary_terms():
    """A BTK-style covalent inhibitor: an acrylamide warhead forms a thioether
    to Cys SG (Michael addition).  The same rule generates the bonded terms
    across the protein-ligand boundary with no per-chemistry code."""
    edges = [
        ("CA", "CB"), ("CB", "SG"), ("CB", "HB1"), ("CB", "HB2"),   # Cys
        ("CL", "OL"), ("CL", "NL"), ("CL", "CY"), ("CY", "HY"),     # ligand amide
        ("CY", "CZ"), ("CZ", "HZ1"), ("CZ", "HZ2"),                 # warhead beta-C
        ("SG", "CZ"),                                               # forming bond
    ]
    atoms, adj = _graph(edges)
    m = induced_bond_morph([(atoms["SG"], atoms["CZ"])], adj, "form")
    bonds = {_names(i) for i in m.interactions if i.kind == "bond"}
    angs = {_canon_ang(_names(i)) for i in m.interactions if i.kind == "angle"}
    dihs = {_canon_dih(_names(i)) for i in m.interactions if i.kind == "dihedral"}

    assert bonds == {("SG", "CZ")}
    assert _canon_ang(("CB", "SG", "CZ")) in angs      # protein-side angle
    assert _canon_ang(("SG", "CZ", "CY")) in angs      # ligand-side angle
    assert _canon_dih(("CA", "CB", "SG", "CZ")) in dihs   # into the protein
    assert _canon_dih(("CL", "CY", "CZ", "SG")) in dihs   # into the ligand
