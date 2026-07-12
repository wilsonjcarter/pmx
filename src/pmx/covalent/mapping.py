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
    from ..ffparser import RTPParser
    from ..utils import get_ff_path
    ffpath = get_ff_path(ff)
    rtp_file = os.path.join(ffpath, "merged.rtp")
    if not os.path.isfile(rtp_file):
        rtp_file = os.path.join(ffpath, "aminoacids.rtp")
    rtp = RTPParser(rtp_file)
    if residue.resname not in rtp:
        raise KeyError(f"{residue.resname} not in {rtp_file}")
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


def _graph_heavy_pairs(res1, res2) -> List[Tuple[int, int]]:
    """Heavy-atom correspondence from a 2-D MCS, disambiguated by atom name."""
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

    # choose the (match1, match2) alignment maximizing atom-name agreement
    def name_score(m1, m2):
        return sum(
            1 for a, b in zip(m1, m2)
            if _canon_name(res1.atoms[heavy1[a]].name)
            == _canon_name(res2.atoms[heavy2[b]].name)
        )

    best, best_score = None, -1
    for m1 in matches1:
        for m2 in matches2:
            s = name_score(m1, m2)
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
    kept = [(i, j) for i, j in heavy_pairs if i not in fi and j not in fj]
    for i, _ in forced:
        method[i] = "backbone"
    return forced + kept


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
    method : {'auto', 'graph', 'name'}
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

    use_graph = (method == "graph") or (method == "auto" and can_graph)

    if use_graph:
        heavy_pairs = _graph_heavy_pairs(res1, res2)
        for i, _ in heavy_pairs:
            method_used[i] = "mcs"
        heavy_pairs = _force_backbone(res1, res2, heavy_pairs, method_used)
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
