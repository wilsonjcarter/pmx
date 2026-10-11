"""Topology-first atom mapping for covalent residue mutations.

The classic pmx ligand mapping (``pmx.ligand_alchemy.LigandAtomMapping``) targets
*arbitrary* ligand pairs with a 3-D maximum-common-substructure search plus
distance-based pruning, which makes the mapping sensitive to embedding/alignment
quality.  (The former ``pmx.nsaa.domapping`` 3-D residue mapper, once used by
``mutate_nsaa``, has been retired in favour of the graph mapping below.)

:func:`map_residues` instead maps on the **molecular graph**, so a pair is only
formed when the two atoms occupy corresponding positions in the connectivity:

* **Graph MCS core** (default) -- a 2-D maximum common substructure over the
  *heavy-atom* graph (element + connectivity, no coordinates).  This is the
  correctness guarantee: matched atoms provably share local topology, so an
  atom that merely *shares a name* but sits elsewhere in the graph cannot be
  mismatched.  Hydrogens are then attached to their mapped heavy parent.
* **Atom names as a prior** -- when the MCS admits several symmetric matches
  (e.g. equivalent branches, interchangeable hydrogens), the one maximizing
  atom-name agreement is chosen, giving a deterministic, chemically-intended
  mapping.  Names never *create* a pair the graph does not support.

Unmatched atoms fall out as dummies: WT-only -> *disappearing* (real in state
A, dummy in B), target-only -> *appearing* (dummy in A, real in B).

A name-only fallback (``method='name'``) is kept for residues without
connectivity, but it is unsafe against name collisions and warns accordingly.

The graph path requires RDKit and that both residues carry bond connectivity
(``atom.bonds``); :func:`attach_bonds_from_rtp` populates it from a force
field's RTP for an arbitrary residue.
"""

from __future__ import annotations

import logging
import os
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    from rdkit.Chem import rdFMCS
    _RDKIT_AVAILABLE = True
except ImportError:  # pragma: no cover
    _RDKIT_AVAILABLE = False


# ── name / element helpers ────────────────────────────────────────────────────

def _canon_name(name: str) -> str:
    """Canonicalize atom names across PDB/GROMACS hydrogen conventions.

    A leading digit (PDB style, ``1HB``) is rotated to the end to match GROMACS
    style (``HB1``).  Non-hydrogen names are unaffected in practice.
    """
    n = name.strip()
    if len(n) > 1 and n[0].isdigit():
        return n[1:] + n[0]
    return n


def _element_of(atom) -> str:
    """Best-effort chemical element for a residue atom."""
    sym = (getattr(atom, "symbol", "") or "").strip()
    if sym:
        return sym[0].upper() + sym[1:].lower()
    n = atom.name.strip()
    if n and n[0].isdigit():
        n = n[1:]
    return n[0].upper() if n else "C"


def _is_h(atom) -> bool:
    return _element_of(atom) == "H"


def _has_bonds(res) -> bool:
    return any(getattr(a, "bonds", None) for a in res.atoms)


# ── result container ──────────────────────────────────────────────────────────

@dataclass
class ResidueMapping:
    """Correspondence between a WT residue (state A) and a target (state B).

    ``pairs`` are ``(i, j)`` index pairs into ``res1.atoms`` / ``res2.atoms``.
    ``disappear`` are res1 indices with no partner (real in A, dummy in B);
    ``appear`` are res2 indices with no partner (dummy in A, real in B).
    ``method`` records how each res1 atom was matched (``'mcs'`` / ``'h-parent'``
    / ``'name'``).
    """

    res1: object
    res2: object
    pairs: List[Tuple[int, int]] = field(default_factory=list)
    disappear: List[int] = field(default_factory=list)
    appear: List[int] = field(default_factory=list)
    method: Dict[int, str] = field(default_factory=dict)

    def pairs1(self) -> List[int]:
        return [i for i, _ in self.pairs]

    def pairs2(self) -> List[int]:
        return [j for _, j in self.pairs]

    def named_pairs(self) -> List[Tuple[str, str]]:
        return [(self.res1.atoms[i].name, self.res2.atoms[j].name)
                for i, j in self.pairs]

    def summary(self) -> str:
        return (f"{len(self.pairs)} mapped, "
                f"{len(self.disappear)} disappear, {len(self.appear)} appear")


# ── connectivity from RTP (for arbitrary residues) ────────────────────────────

def attach_bonds_from_rtp(residue, ff: str) -> None:
    """Populate ``atom.bonds`` for *residue* from the force field's RTP entry.

    Use this to give a WT residue (e.g. loaded from a protein) the connectivity
    the graph mapper needs.  Intra-residue bonds only; ``+``/``-`` inter-residue
    bonds are skipped.  No-op for atoms already bonded.
    """
    from ..ffparser import load_ff_rtp
    from ..utils import get_ff_path
    ffpath = get_ff_path(ff)
    rtp = load_ff_rtp(ff)      # base RTP + any topic RTPs
    if residue.resname not in rtp:
        searched = ", ".join(os.path.basename(s) for s in rtp.sources)
        raise KeyError(f"{residue.resname} not in {ffpath} (searched: {searched})")
    by_name = {}
    for a in residue.atoms:
        by_name[a.name] = a
        by_name[_canon_name(a.name)] = a
    for a1, a2 in rtp[residue.resname]["bonds"]:
        if a1.startswith(("+", "-")) or a2.startswith(("+", "-")):
            continue
        x = by_name.get(a1) or by_name.get(_canon_name(a1))
        y = by_name.get(a2) or by_name.get(_canon_name(a2))
        if x is not None and y is not None and y not in x.bonds:
            x.bonds.append(y)
            y.bonds.append(x)


# ── graph MCS ─────────────────────────────────────────────────────────────────

def _heavy_mol(res):
    """Build an RDKit mol of the heavy atoms + heavy-heavy single bonds.

    Returns (mol, heavy_indices) where mol atom k corresponds to
    res.atoms[heavy_indices[k]].
    """
    heavy = [i for i, a in enumerate(res.atoms) if not _is_h(a)]
    pos = {ri: k for k, ri in enumerate(heavy)}
    rw = Chem.RWMol()
    for ri in heavy:
        rw.AddAtom(Chem.Atom(_element_of(res.atoms[ri])))
    seen = set()
    for ri in heavy:
        for nb in getattr(res.atoms[ri], "bonds", []):
            rj = _index_of(res, nb)
            if rj is None or rj not in pos:
                continue
            key = tuple(sorted((ri, rj)))
            if key in seen:
                continue
            seen.add(key)
            rw.AddBond(pos[ri], pos[rj], Chem.BondType.SINGLE)
    mol = rw.GetMol()
    Chem.SanitizeMol(
        mol,
        sanitizeOps=(Chem.SanitizeFlags.SANITIZE_ALL
                     ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES),
    )
    return mol, heavy


def _index_of(res, atom) -> Optional[int]:
    for i, a in enumerate(res.atoms):
        if a is atom:
            return i
    return None


def _heavy_neighbors(res, i) -> List[int]:
    """Indices of *i*'s heavy-atom neighbours in *res*."""
    out = []
    for b in getattr(res.atoms[i], "bonds", []):
        k = _index_of(res, b)
        if k is not None and not _is_h(res.atoms[k]):
            out.append(k)
    return out


def _ring_size(res, i) -> int:
    """Size of the SMALLEST cycle through heavy atom *i*; 0 if acyclic.

    Membership alone is too weak a constraint: PHE's CD1/CE1/CZ sit in a
    six-ring and TRP's CD1/NE1/CE2 in the fused five-ring, so an "is it in a
    ring" test let them pair -- including a carbon onto a nitrogen -- which would
    put both a 5-cycle and a 6-cycle in the hybrid's single bond list and make
    ring connectivity state-dependent.
    """
    nb = _heavy_neighbors(res, i)
    best = 0
    for a in range(len(nb) - 1):
        for b in range(a + 1, len(nb)):
            # shortest path nb[a] -> nb[b] avoiding i, by BFS
            start, goal = nb[a], nb[b]
            seen, frontier, dist = {i, start}, [start], 1
            while frontier and not best or frontier:
                nxt = []
                for cur in frontier:
                    if cur == goal:
                        cyc = dist + 1
                        best = cyc if not best else min(best, cyc)
                        nxt = []
                        break
                    for k in _heavy_neighbors(res, cur):
                        if k != i and k not in seen:
                            seen.add(k)
                            nxt.append(k)
                else:
                    frontier, dist = nxt, dist + 1
                    continue
                break
    return best


# ── prototype: anchored shell matcher ─────────────────────────────────────────
# An alternative to FindMCS for the RESIDUE case, where a backbone anchor exists.
# Our problem is not maximum-common-substructure: the backbone is forced, a pair
# must sit at equal bond depth from CA, and a pair's parents must correspond, so
# the correspondence is determined shell by shell outward from CA and each shell
# has a handful of candidates (TRP, the largest side chain, has 10 heavy atoms).
#
# FindMCS optimises substructure SIZE, which is the wrong objective here: for
# PHE->TRP it matched PHE's six-ring onto TRP's benzo ring (6 atoms, a legal
# complete-ring match, beating the correct CB+CG=2) and _enforce_anchor_depth then
# had to discard the lot, taking CB with it.  It also cannot express the preference
# we want -- its element handling is binary (CompareElements forbids a
# cross-element map, CompareAny is indifferent) -- whereas upstream pmx maps
# C<->O, C<->S, O<->S and N<->O when no same-element candidate exists at that
# position.  A cost says "prefer same element, allow different" in one line.
#
# Every term is swap-invariant, so the matching is identical in both directions.
_COST_CROSS_ELEMENT = 10.0
# Carbon-onto-heteroatom is a bigger chemical change than heteroatom-onto-
# heteroatom, and without separating them the two tie: THR has OG1 and CG2 at
# depth 2 against CYS's SG, and a flat cross-element cost let CG2->SG (methyl
# carbon onto thiol sulfur) win on the name/degree tiebreak.  Upstream maps
# OG1->SG -- hydroxyl oxygen onto thiol sulfur, both bearing their H.
_COST_CARBON_FOR_HETERO = 5.0
_COST_NAME_MISMATCH = 1.0
_COST_DEGREE = 0.1
# The ring ATTACHMENT atom may map across a ring-size change, which is what
# upstream does (PHE CG -> TRP CG, CA -> CY in F2W).  Measured: it resolves
# exactly the 16 aromatic pairs, 1 atom and 2 dummies cheaper each, with the
# three mapping invariants still holding on all 380 ordered pairs.
_attach_across_rings = True


def _is_ring_attachment(res, i, mapped, ringsz) -> bool:
    """*i* is a ring atom whose mapped neighbours are all ACYCLIC."""
    nb = [p for p in _heavy_neighbors(res, i) if p in mapped]
    return bool(nb) and all(not ringsz.get(p) for p in nb)


def _has_anchor(res) -> bool:
    """Does *res* carry the CA the rooted matcher roots on?"""
    return any(_canon_name(a.name) == _ANCHOR for a in res.atoms)


def _anchored_pairs(res1, res2, method):
    """Rooted shell-by-shell heavy-atom correspondence, CA outward."""
    d1, d2 = _graph_depth(res1), _graph_depth(res2)
    if not d1 or not d2:
        return []

    pairs = _force_backbone(res1, res2, [], method)      # anchor, by name
    mapped = dict(pairs)
    used2 = set(mapped.values())

    heavy1 = [i for i, a in enumerate(res1.atoms) if not _is_h(a)]
    heavy2 = [j for j, a in enumerate(res2.atoms) if not _is_h(a)]
    ring1 = {i: _ring_size(res1, i) for i in heavy1}
    ring2 = {j: _ring_size(res2, j) for j in heavy2}

    depths = sorted({d for d in d1.values() if d is not None and d > 0})
    for d in depths:
        while True:
            rev = {j: i for i, j in mapped.items()}
            c1 = [i for i in heavy1 if i not in mapped and d1.get(i) == d]
            c2 = [j for j in heavy2 if j not in used2 and d2.get(j) == d]
            if not c1 or not c2:
                break

            def cost(i, j):
                # HARD: smallest-ring size must agree (0 == acyclic).  A hybrid
                # carries one bond list, so pairing a five-ring atom with a
                # six-ring atom would make ring connectivity state-dependent.
                if ring1[i] != ring2[j]:
                    # Exception, opt-in: the ring ATTACHMENT atom, whose only
                    # mapped neighbour is acyclic (CB).  Upstream maps PHE CG
                    # onto TRP CG across a six-/five-ring change; the bond that
                    # persists is CG-CB, which is outside both rings.
                    if not (_attach_across_rings
                            and ring1[i] and ring2[j]
                            and _is_ring_attachment(res1, i, mapped, ring1)
                            and _is_ring_attachment(res2, j,
                                                    {v: k for k, v in mapped.items()},
                                                    ring2)):
                        return None
                # HARD: the parents must already correspond, which is what keeps
                # the mapping connected to CA and the dummies terminal.
                p1 = [p for p in _heavy_neighbors(res1, i) if p in mapped]
                if not any(mapped[p] in _heavy_neighbors(res2, j) for p in p1):
                    return None
                n1 = _canon_name(res1.atoms[i].name)
                n2 = _canon_name(res2.atoms[j].name)
                e1 = _element_of(res1.atoms[i])
                e2 = _element_of(res2.atoms[j])
                c = 0.0
                if e1 != e2:
                    c += _COST_CROSS_ELEMENT
                    if ('C' in (e1, e2)) and not {e1, e2} <= {'C'}:
                        c += _COST_CARBON_FOR_HETERO
                if n1 != n2:
                    c += _COST_NAME_MISMATCH
                c += _COST_DEGREE * abs(len(_heavy_neighbors(res1, i))
                                        - len(_heavy_neighbors(res2, j)))
                return (c, tuple(sorted((n1, n2))))

            opts = {}
            for i in c1:
                for j in c2:
                    k = cost(i, j)
                    if k is not None:
                        opts[(i, j)] = k
            if not opts:
                break
            best1, best2 = {}, {}
            for (i, j), k in sorted(opts.items(), key=lambda t: (t[1], t[0])):
                best1.setdefault(i, (k, j))
                best2.setdefault(j, (k, i))
            took = False
            for i, (k, j) in sorted(best1.items()):
                if best2.get(j, (None, None))[1] != i:
                    continue                      # not mutually best -> symmetric
                if i in mapped or j in used2:
                    continue
                mapped[i] = j
                used2.add(j)
                method[i] = ("isosteric" if _element_of(res1.atoms[i])
                             != _element_of(res2.atoms[j]) else "anchored")
                took = True
            if not took:
                break
    return sorted(mapped.items())


def _graph_heavy_pairs(res1, res2) -> List[Tuple[int, int]]:
    """Heavy-atom correspondence from a 2-D MCS, disambiguated by position."""
    mol1, heavy1 = _heavy_mol(res1)
    mol2, heavy2 = _heavy_mol(res2)
    if mol1.GetNumAtoms() == 0 or mol2.GetNumAtoms() == 0:
        return []

    res = rdFMCS.FindMCS(
        [mol1, mol2],
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareAny,
        ringMatchesRingOnly=True,
        completeRingsOnly=True,
        timeout=10,
    )
    if not res.smartsString:
        return []
    q = Chem.MolFromSmarts(res.smartsString)
    if q is None:
        return []
    matches1 = mol1.GetSubstructMatches(q, uniquify=False, maxMatches=2000)
    matches2 = mol2.GetSubstructMatches(q, uniquify=False, maxMatches=2000)
    if not matches1 or not matches2:
        return []

    # Choose the (match1, match2) alignment by DEPTH agreement first, name second.
    # Depth from CA is the positional prior and must outrank the label: with the
    # backbone excluded from the MCS there is no N-CA anchor left to pin the side
    # chain, so a pure name score picks a target atom that merely reuses a name at
    # a different graph position (WT N-CA-CB vs target N-CA-CG-CB would align
    # CB<->CB rather than the positionally correct CB<->CG).  Scoring depth first
    # restores the "map by position, not by name" guarantee, and name agreement
    # still breaks ties among equal-depth alignments (symmetric substituents).
    dd1, dd2 = _graph_depth(res1), _graph_depth(res2)

    def align_score(m1, m2):
        depth = name = 0
        for a, b in zip(m1, m2):
            i, j = heavy1[a], heavy2[b]
            if dd1.get(i) is not None and dd1.get(i) == dd2.get(j):
                depth += 1
            if (_canon_name(res1.atoms[i].name)
                    == _canon_name(res2.atoms[j].name)):
                name += 1
        return (depth, name)

    best, best_score = None, (-1, -1)
    for m1 in matches1:
        for m2 in matches2:
            s = align_score(m1, m2)
            if s > best_score:
                best_score, best = s, (m1, m2)
    m1, m2 = best
    return [(heavy1[a], heavy2[b]) for a, b in zip(m1, m2)]


# The peptide backbone is the invariant scaffold of any residue mutation, so
# these always correspond by name -- even when one side puts them in a ring
# (proline), where a pure ring-matching MCS would otherwise refuse to map them.
_BACKBONE = {"N", "CA", "C", "O", "OT1", "OT2", "OC1", "OC2", "OXT"}


def _force_backbone(res1, res2, heavy_pairs, method):
    """Force backbone atoms to map by name, overriding the MCS where needed."""
    n1 = {_canon_name(a.name): i for i, a in enumerate(res1.atoms)}
    n2 = {_canon_name(a.name): j for j, a in enumerate(res2.atoms)}
    forced = []
    for nm in _BACKBONE:
        if nm in n1 and nm in n2:
            forced.append((n1[nm], n2[nm]))
    fi = {i for i, _ in forced}
    fj = {j for _, j in forced}
    # A backbone atom must never map onto a NON-backbone one, and vice versa.
    # _force_backbone can only pin a name present in BOTH residues, so a target
    # supplied as a capped side-chain fragment (SFR: CA, CB, OG + a CH3 stand-in,
    # no N/C/O) leaves the WT's C and O free, and the MCS slid the backbone onto
    # the fragment -- C onto CB, O onto OG.  Depth cannot see it either, because
    # the cap makes the two branches off CA symmetric.  One rule, no new machinery.
    def _cross_class(i, j):
        return ((_canon_name(res1.atoms[i].name) in _BACKBONE)
                != (_canon_name(res2.atoms[j].name) in _BACKBONE))
    kept = [(i, j) for i, j in heavy_pairs
            if i not in fi and j not in fj and not _cross_class(i, j)]
    for i, _ in forced:
        method[i] = "backbone"
    return forced + kept


# The alchemical anchor.  Forcing the backbone to map by name (above) pins CA, but it
# does NOT constrain what the atoms hanging off CA map to: FindMCS maximises
# common-substructure SIZE, and _force_backbone only discards pairs that collide with
# a forced index.  So the MCS is free to match from the TIP of a side chain and slide
# the whole chain by a bond -- ATB(CB) pairs with HFL(CG) because that lines the two
# CF3 groups up and buys more matched atoms than the chemically correct CB<->CB.
_ANCHOR = "CA"


def _graph_depth(res, root: str = _ANCHOR) -> Dict[int, int]:
    """Bond-graph distance of every atom from *root*; ``{}`` if *root* is absent."""
    idx = {}
    for i, a in enumerate(res.atoms):
        idx.setdefault(_canon_name(a.name), i)
    if root not in idx:
        return {}
    adj: Dict[int, set] = {i: set() for i in range(len(res.atoms))}
    for i, a in enumerate(res.atoms):
        for nb in (getattr(a, "bonds", None) or []):
            j = _index_of(res, nb)
            if j is not None:
                adj[i].add(j)
                adj[j].add(i)
    depth = {idx[root]: 0}
    queue = [idx[root]]
    while queue:
        u = queue.pop(0)
        for v in adj[u]:
            if v not in depth:
                depth[v] = depth[u] + 1
                queue.append(v)
    return depth


def _enforce_anchor_depth(res1, res2, heavy_pairs, method):
    """Drop mapped pairs whose atoms sit at different bond distances from CA.

    A hybrid residue carries ONE bond list, and bond *existence* cannot be made
    state-dependent (only bonded parameters morph).  So a pair at unequal depth is
    un-representable: the atom would have to be bonded straight to CA in one state and
    one bond further out in the other.  The bond-union then resolves the impossible
    constraint differently depending on which residue is state A -- an extra bond on CA
    (five substituents) when the shorter residue is A, or a side chain whose only path
    to the backbone runs through a dummy when the longer one is.  Either way the TARGET
    end state stops being the target residue, while the source state stays correct.

    Dropping the pair is the safe override: both atoms simply become dummies on their
    own side, which is a valid (if larger) alchemical path.
    """
    d1 = _graph_depth(res1)
    d2 = _graph_depth(res2)
    if not d1 or not d2:
        return heavy_pairs

    kept, dropped = [], []
    for i, j in heavy_pairs:
        if method.get(i) == "backbone":
            kept.append((i, j))
            continue
        a, b = d1.get(i), d2.get(j)
        if a is None or b is None or a == b:
            kept.append((i, j))
        else:
            dropped.append((i, j, a, b))

    # Dropping a pair can orphan the ones beyond it: PHE->TRP keeps the slid ring
    # tip (CZ/CE2/CD2) after its root is discarded, leaving a mapped island whose
    # only route to CA runs through a dummy.  A dummy carries full bonded terms, so
    # the molecule is not severed, but a real fragment tethered by a ghost is not a
    # chemically sensible end state and it breaks the "dummies are terminal" rule.
    # Sweep them, repeatedly, so the mapped set is always connected to CA.
    ca = next((i for i, a in enumerate(res1.atoms)
               if _canon_name(a.name) == _ANCHOR), None)
    if ca is not None:
        while True:
            mp = {i for i, _ in kept}
            if ca not in mp:
                break
            seen, stack = {ca}, [ca]
            while stack:
                for k in _heavy_neighbors(res1, stack.pop()):
                    if k in mp and k not in seen:
                        seen.add(k)
                        stack.append(k)
            orphans = {i for i in mp
                       if i not in seen and not _is_h(res1.atoms[i])}
            if not orphans:
                break
            for i, j in [p for p in kept if p[0] in orphans]:
                logger.warning("map_residues:   orphan %s dropped (reaches %s "
                               "only through a dummy)", res1.atoms[i].name, _ANCHOR)
                dropped.append((i, j, -1, -1))
                method.pop(i, None)
            kept = [p for p in kept if p[0] not in orphans]

    if dropped:
        logger.warning(
            "map_residues: MCS anchor break detected, %s -> %s: %d of %d mapped "
            "pairs sit at different bond distances from %s, so the MCS has matched "
            "from the tip and slid the chain. Overriding: these pairs are dropped and "
            "the atoms become dummies on their own side. Keeping them would leave the "
            "target end state not equal to the target residue, because a hybrid has a "
            "single bond list and bond existence is not state-dependent.",
            getattr(res1, "resname", "?"), getattr(res2, "resname", "?"),
            len(dropped), len(heavy_pairs), _ANCHOR)
        for i, j, a, b in dropped:
            logger.warning("map_residues:   override %s (d=%d) <-> %s (d=%d)",
                           res1.atoms[i].name, a, res2.atoms[j].name, b)
            method.pop(i, None)
    return kept


def _h_neighbors(res, i) -> List[int]:
    out = []
    for nb in getattr(res.atoms[i], "bonds", []):
        j = _index_of(res, nb)
        if j is not None and _is_h(res.atoms[j]):
            out.append(j)
    return out


def _map_hydrogens(res1, res2, heavy_pairs, method):
    """Attach hydrogens to their mapped heavy parent, paired by canonical name."""
    h_pairs = []
    for i, j in heavy_pairs:
        h1 = _h_neighbors(res1, i)
        h2 = _h_neighbors(res2, j)
        used = set()
        # first by canonical name, then greedily by remaining order
        for a in list(h1):
            ca = _canon_name(res1.atoms[a].name)
            match = next((b for b in h2 if b not in used
                          and _canon_name(res2.atoms[b].name) == ca), None)
            if match is not None:
                h_pairs.append((a, match))
                used.add(match)
                method[a] = "h-parent"
                h1.remove(a)
        for a, b in zip(h1, [b for b in h2 if b not in used]):
            h_pairs.append((a, b))
            used.add(b)
            method[a] = "h-parent"
    return h_pairs


# ── name-only fallback ────────────────────────────────────────────────────────

def _name_pairs(res1, res2, element_guard, method):
    r2_by_canon: Dict[str, List[int]] = {}
    for j, a in enumerate(res2.atoms):
        r2_by_canon.setdefault(_canon_name(a.name), []).append(j)
    used2: set = set()
    pairs: List[Tuple[int, int]] = []
    for i, a in enumerate(res1.atoms):
        canon = _canon_name(a.name)
        chosen = None
        for j in r2_by_canon.get(canon, []):
            if j in used2:
                continue
            if element_guard and _element_of(a) != _element_of(res2.atoms[j]):
                continue
            chosen = j
            break
        if chosen is not None:
            pairs.append((i, chosen))
            used2.add(chosen)
            method[i] = "name"
    return pairs


# ── public API ────────────────────────────────────────────────────────────────

def map_residues(res1, res2, method: str = "auto",
                 element_guard: bool = True) -> ResidueMapping:
    """Map WT residue *res1* (state A) onto target *res2* (state B).

    Parameters
    ----------
    res1, res2 : residue-like
        Objects with an ``.atoms`` list of pmx :class:`~pmx.atom.Atom`.  For the
        graph method both must carry connectivity (``atom.bonds``); see
        :func:`attach_bonds_from_rtp`.
    method : {'auto', 'anchored', 'graph', 'name'}
        ``'anchored'`` -- rooted shell matcher from CA (the default via ``'auto'``).
        ``'graph'`` -- MCS on the molecular graph (safe; requires RDKit + bonds).
        ``'name'``  -- name+element anchor only (unsafe against name collisions).
        ``'auto'``  -- graph when possible, else name with a warning.
    element_guard : bool
        For the name path, require matched names to share an element.
    Returns
    -------
    ResidueMapping
    """
    method_used: Dict[int, str] = {}
    can_graph = _RDKIT_AVAILABLE and _has_bonds(res1) and _has_bonds(res2)

    if method == "graph" and not can_graph:
        raise ValueError(
            "method='graph' requires RDKit and bond connectivity on both "
            "residues; attach bonds (e.g. attach_bonds_from_rtp) or use "
            "method='name'.")

    # 'auto' prefers the rooted shell matcher.  It beats FindMCS on every measure
    # we have: agreement with the shipped library 198 -> 16 differing hybrids of
    # 650, identical in both force fields, with end states 380/380, reverse
    # symmetry 380/380 and the three mapping invariants clean.  FindMCS optimises
    # substructure SIZE, which is the wrong objective once the backbone is forced
    # and depth is constrained, and it cannot express "prefer the same element,
    # allow a different one" at all.  Pass method='graph' for FindMCS explicitly.
    if method == "auto" and can_graph:
        method = "anchored"

    if method == "anchored" and not (_has_anchor(res1) and _has_anchor(res2)):
        # No CA on one side, so there is no root and no depth coordinate: a GFP
        # chromophore (CRO has CA1/CA2/CA3 and no CA, because it replaces the
        # backbone) or an arbitrary ligand.  The rooted matcher would decouple
        # everything, and maximum-common-substructure is genuinely the right
        # primitive for an UNANCHORED pair, so hand those back to FindMCS.
        logger.info("map_residues: no %s anchor on both residues (%s -> %s); "
                    "using the MCS path", _ANCHOR,
                    getattr(res1, "resname", "?"), getattr(res2, "resname", "?"))
        method = "auto"

    if method == "anchored":
        pairs = _anchored_pairs(res1, res2, method_used)
        pairs = pairs + _map_hydrogens(res1, res2, pairs, method_used)
        matched1 = {i for i, _ in pairs}
        matched2 = {j for _, j in pairs}
        return ResidueMapping(
            res1=res1, res2=res2, pairs=pairs,
            disappear=[i for i in range(len(res1.atoms)) if i not in matched1],
            appear=[j for j in range(len(res2.atoms)) if j not in matched2],
            method=method_used)

    use_graph = (method == "graph") or (method == "auto" and can_graph)

    if use_graph:
        heavy_pairs = _graph_heavy_pairs(res1, res2)
        for i, _ in heavy_pairs:
            method_used[i] = "mcs"
        heavy_pairs = _force_backbone(res1, res2, heavy_pairs, method_used)
        heavy_pairs = _enforce_anchor_depth(res1, res2, heavy_pairs, method_used)
        h_pairs = _map_hydrogens(res1, res2, heavy_pairs, method_used)
        pairs = heavy_pairs + h_pairs
    else:
        if method == "auto":
            logger.warning(
                "map_residues: no connectivity/RDKit available; falling back "
                "to name-only matching, which is unsafe against atom-name "
                "collisions. Attach bonds for a graph-verified mapping.")
        pairs = _name_pairs(res1, res2, element_guard, method_used)

    matched1 = {i for i, _ in pairs}
    matched2 = {j for _, j in pairs}
    disappear = [i for i in range(len(res1.atoms)) if i not in matched1]
    appear = [j for j in range(len(res2.atoms)) if j not in matched2]

    return ResidueMapping(res1=res1, res2=res2, pairs=pairs,
                          disappear=disappear, appear=appear,
                          method=method_used)
