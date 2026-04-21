"""Non-standard amino acid (NSAA / PTM) mutation utilities.

Two-step workflow
-----------------
1. **Structure prep** — call :func:`mutate_nsaa` to build the hybrid residue
   PDB, RTP, and MTP files needed by pdb2gmx.

2. **Topology filling** — call :func:`~pmx.alchemy.gen_hybrid_top` with the
   ``extra_mtp_files`` argument pointing at the MTP generated in step 1.

Typical usage::

    from pmx.model import Model
    from pmx.nsaa import mutate_nsaa
    from pmx.alchemy import gen_hybrid_top
    from pmx.forcefield import Topology

    m = Model("protein.pdb", bPDBTER=True, rename_atoms=True, scale_coords='A')
    mutate_nsaa(m, mut_resid=42, mut_resname="SEP", ff="charmm36m-mut",
                nonstandardPDB="sep.pdb", itp="sep.itp", inplace=True)
    m.write("mutant.pdb")

    # After pdb2gmx:
    top = Topology("topol.top", ff="charmm36m-mut", version="new")
    pmxtop, itps = gen_hybrid_top(top, extra_mtp_files=["A2SEP.mtp"])
"""

from __future__ import annotations

import logging
import os
import re
import sys
from copy import deepcopy
from typing import Dict, List, Optional, Tuple

import numpy as np

from . import library
from .atom import Atom
from .geometry import Rotation
from .model import Model
from .parser import kickOutComments
from .utils import get_ff_path, list2file
from .library import pmx_data_file

# RDKit is an optional dependency — only needed for domapping()
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdmolops
    from rdkit.Chem.AllChem import AlignMol, EmbedMultipleConfs
    from .ligand_alchemy import LigandAtomMapping
    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False

logger = logging.getLogger(__name__)

# ── Module-level constants ────────────────────────────────────────────────────

# Inverted alias table: RTP atom name → PDB atom name (per residue)
_inverted_aliases: Dict[str, Dict[str, str]] = {}
for _k in library._aliases:
    _inverted_aliases[_k] = {v: k for k, v in library._aliases[_k].items()}
_inverted_aliases.setdefault("ILE", {})["CD1"] = "CD1"

# Standard atomic masses for element-symbol fallback
_ELEM_MASS = {
    'H': 1.008,  'C': 12.011, 'N': 14.007, 'O': 15.999,
    'S': 32.06,  'P': 30.974, 'F': 18.998, 'CL': 35.45,
    'BR': 79.904,'I': 126.904,'SE': 78.971,
}

# Backbone atom names (used for parameter assignment and backbone cloning)
_BACKBONE = {'N', 'HN', 'H', 'CA', 'HA', 'C', 'O', 'OT1', 'OT2'}
_N_SIDE   = {'N', 'HN', 'H', 'HT1', 'HT2', 'HT3', 'H1', 'H2', 'H3'}
_C_SIDE   = {'C', 'O', 'OT1', 'OT2'}

_CHARMM_RTP_HEADER = """\
[ bondedtypes ]
; Col 1: Type of bond
; Col 2: Type of angles
; Col 3: Type of proper dihedrals
; Col 4: Type of improper dihedrals
; Col 5: Generate all dihedrals if 1, only heavy atoms of 0.
; Col 6: Number of excluded neighbors for nonbonded interactions
; Col 7: Generate 1,4 interactions between pairs of hydrogens if 1
; Col 8: Remove propers over the same bond as an improper if it is 1
; bonds  angles  dihedrals  impropers  all_dihedrals  nrexcl  HH14  RemoveDih
    1       5        9          2            1           3      1       0
"""

_AMBER_RTP_HEADER = """\
[ bondedtypes ]
; Col 1: Type of bond
; Col 2: Type of angles
; Col 3: Type of proper dihedrals
; Col 4: Type of improper dihedrals
; Col 5: Generate all dihedrals if 1, only heavy atoms of 0.
; Col 6: Number of excluded neighbors for nonbonded interactions
; Col 7: Generate 1,4 interactions between pairs of hydrogens if 1
; Col 8: Remove impropers over the same bond as a proper if it is 1
; bonds  angles  dihedrals  impropers all_dihedrals nrexcl HH14 RemoveDih
     1       1          9          4        1         3      1     0
"""

# CHARMM atom-renaming specification for pdb2gmx
_ARN_CONTENT = """\
; atom renaming specification
; residue gromacs    forcefield
*          H          HN
"""


# =============================================================================
# ITP parsing
# =============================================================================

def _parse_itp_atoms(itp_file: str) -> List[dict]:
    """Parse the ``[ atoms ]`` section from a GROMACS ITP file.

    Returns a list of dicts with keys: ``nr``, ``type``, ``name``,
    ``charge``, ``mass``.
    """
    atoms = []
    in_atoms = False
    with open(itp_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(';'):
                continue
            if line.startswith('['):
                in_atoms = (line.strip('[]').strip() == 'atoms')
                continue
            if in_atoms:
                parts = line.split()
                if len(parts) >= 7:
                    try:
                        atoms.append({
                            'nr':     int(parts[0]),
                            'type':   parts[1],
                            'name':   parts[4],
                            'charge': float(parts[6]),
                            'mass':   float(parts[7]) if len(parts) > 7 else 0.0,
                        })
                    except (ValueError, IndexError):
                        continue
    return atoms


def _parse_itp_bonds(itp_file: str) -> List[Tuple[int, int]]:
    """Parse the ``[ bonds ]`` section from a GROMACS ITP file.

    Returns a list of ``(ai, aj)`` tuples of 1-based atom numbers.
    """
    bonds = []
    in_bonds = False
    with open(itp_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(';'):
                continue
            if line.startswith('['):
                in_bonds = (line.strip('[]').strip() == 'bonds')
                continue
            if in_bonds:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        bonds.append((int(parts[0]), int(parts[1])))
                    except ValueError:
                        continue
    return bonds


def _parse_itp_section(itp_file: str, section_name: str,
                       n_index_cols: int) -> List[dict]:
    """Generic parser for numbered ITP sections (bonds, angles, dihedrals).

    Returns a list of dicts::

        {'indices': [int, ...], 'funct': int, 'params': [str, ...]}

    Inline comments (semicolons) are stripped before parsing.
    """
    result = []
    in_section = False
    with open(itp_file) as fh:
        for raw in fh:
            line = raw.split(';')[0].strip()
            if not line:
                continue
            if line.startswith('['):
                in_section = (line.strip('[]').strip() == section_name)
                continue
            if not in_section:
                continue
            parts = line.split()
            if len(parts) < n_index_cols + 1:
                continue
            try:
                indices = [int(parts[i]) for i in range(n_index_cols)]
                funct   = int(parts[n_index_cols])
            except ValueError:
                continue
            result.append({
                'indices': indices,
                'funct':   funct,
                'params':  parts[n_index_cols + 1:],
            })
    return result


def _parse_itp_bonds_full(itp_file: str) -> List[dict]:
    """Parse ``[ bonds ]`` from an ITP file, including explicit parameters."""
    return _parse_itp_section(itp_file, 'bonds', 2)


def _parse_itp_angles(itp_file: str) -> List[dict]:
    """Parse ``[ angles ]`` from an ITP file."""
    return _parse_itp_section(itp_file, 'angles', 3)


def _parse_itp_dihedrals(itp_file: str) -> List[dict]:
    """Parse ``[ dihedrals ]`` from an ITP file (proper and improper)."""
    return _parse_itp_section(itp_file, 'dihedrals', 4)


def _parse_atomtypes_section(atomtypes_itp_file: str) -> List[List[str]]:
    """Parse a ``[ atomtypes ]`` section from a GROMACS ITP file.

    Returns a list of token lists::

        [name, at_num, mass, charge, ptype, sigma, epsilon]

    Only the first ``[ atomtypes ]`` block is read.
    """
    entries = []
    in_section = False
    with open(atomtypes_itp_file) as fh:
        for raw in fh:
            line = raw.split(';')[0].strip()
            if not line:
                continue
            if line.startswith('['):
                section = line.strip('[]').strip()
                in_section = (section == 'atomtypes')
                if in_section:
                    continue
                if entries:   # first block done → stop
                    break
                continue
            if not in_section:
                continue
            parts = line.split()
            if len(parts) >= 7:
                entries.append(parts[:7])
            elif len(parts) >= 6:
                # Some exporters omit the at.num column; pad with '0'
                entries.append([parts[0], '0'] + parts[1:6])
    return entries


def _parse_ffbonded_section(ff_file: str, section_name: str,
                             n_type_cols: int) -> List[dict]:
    """Parse a *types section from ``ffbonded.itp`` (type-name based).

    Returns a list of dicts::

        {'types': [str, ...], 'funct': int, 'params': [str, ...]}
    """
    result = []
    in_section = False
    with open(ff_file) as fh:
        for raw in fh:
            line = raw.split(';')[0].strip()
            if not line:
                continue
            if line.startswith('['):
                in_section = (line.strip('[]').strip() == section_name)
                continue
            if not in_section:
                continue
            parts = line.split()
            if len(parts) < n_type_cols + 1:
                continue
            types = parts[:n_type_cols]
            try:
                funct = int(parts[n_type_cols])
            except ValueError:
                continue
            result.append({
                'types':  types,
                'funct':  funct,
                'params': parts[n_type_cols + 1:],
            })
    return result


def _make_itp_names_unique(itp_atoms: List[dict]) -> Optional[dict]:
    """Rename ITP atom names in-place when they are not all unique.

    Uses *element_symbol + nr* (e.g. ``'C'`` at nr 5 → ``'C5'``).
    This is a no-op when names are already unique.

    Returns ``{nr: unique_name}`` if renaming was done, ``None`` otherwise.
    """
    names = [a['name'] for a in itp_atoms]
    if len(set(names)) == len(names):
        return None
    nr_to_unique = {a['nr']: a['name'] + str(a['nr']) for a in itp_atoms}
    for a in itp_atoms:
        a['name'] = nr_to_unique[a['nr']]
    logger.info("_make_itp_names_unique: renamed %d atoms to unique names.",
                len(itp_atoms))
    return nr_to_unique


# =============================================================================
# Fragment preparation
# =============================================================================

def _find_methyl_cap(itp_atoms: List[dict], itp_bonds: List[Tuple[int, int]]):
    """Identify the CH3 cap group that represents the CA position in an ITP.

    The cap is a carbon with exactly 3 H neighbours and exactly 1 non-H
    neighbour.  When multiple candidates exist the correct one is the CH3
    whose non-H neighbour has the most H atoms — because CA in a capped
    fragment retains extra Hs in place of the missing backbone bonds.

    Returns
    -------
    cap_C_name, cap_H_names, ca_atom_name, ca_H_names
        All ``None`` / ``[]`` / ``None`` / ``[]`` when not found.
    """
    adj = {a['nr']: [] for a in itp_atoms}
    for ai, aj in itp_bonds:
        if ai in adj:
            adj[ai].append(aj)
        if aj in adj:
            adj[aj].append(ai)

    nr_to = {a['nr']: a for a in itp_atoms}

    def _is_h(nr):
        name = nr_to[nr]['name']
        return name.startswith('H') or (
            len(name) > 1 and name[0].isdigit() and name[1].upper() == 'H')

    candidates = []
    for atom in itp_atoms:
        if atom['type'][0].upper() != 'C':
            continue
        nbrs   = adj.get(atom['nr'], [])
        h_nbrs = [n for n in nbrs if _is_h(n)]
        non_h  = [n for n in nbrs if not _is_h(n)]
        if len(h_nbrs) == 3 and len(non_h) == 1:
            ca_nr    = non_h[0]
            ca_h_cnt = sum(1 for n in adj.get(ca_nr, []) if _is_h(n))
            candidates.append((atom, h_nbrs, ca_nr, ca_h_cnt))

    if not candidates:
        logger.warning("Could not auto-detect CH3 cap — no terminal CH3 found.")
        return None, [], None, []

    best = max(candidates, key=lambda x: x[3])
    atom, h_nbrs, ca_nr, _ = best

    if len(candidates) > 1:
        logger.info(
            "Multiple CH3 caps found; selected cap='%s' (neighbour '%s' has "
            "most H atoms).", atom['name'], nr_to[ca_nr]['name'])

    cap_H_names  = [nr_to[n]['name'] for n in h_nbrs]
    ca_atom_name = nr_to[ca_nr]['name']
    ca_H_names   = [nr_to[n]['name'] for n in adj.get(ca_nr, []) if _is_h(n)]
    logger.info("CH3 cap: C='%s', cap_H=%s, CA='%s', CA_H=%s",
                atom['name'], cap_H_names, ca_atom_name, ca_H_names)
    return atom['name'], cap_H_names, ca_atom_name, ca_H_names


def _prepare_nsaa_fragment(residue2, residue1, cap_C_name, cap_H_names,
                           ca_atom_name, ca_H_names=None, term=None):
    """Convert a CH3-capped sidechain fragment into a full residue in-place.

    Steps:

    1. Remove the CH3 cap (``cap_C_name`` + ``cap_H_names``) from *residue2*.
    2. Rename ``ca_atom_name`` → ``'CA'``; keep one CA hydrogen as ``'HA'``;
       drop any additional CA hydrogens (the capped fragment has one extra H
       where the backbone bonds would normally be).
    3. Clone backbone atoms from *residue1* via SVD rigid-body superposition
       anchored on CB (preferred) or CA+HA.

    The final atom order is:
    ``[N-side backbone] + [CA + sidechain] + [C-side backbone]``.
    """
    # ── Step 1: remove cap ────────────────────────────────────────────────────
    cap_names = set(cap_H_names)
    if cap_C_name:
        cap_names.add(cap_C_name)
    remove_ids = {id(a) for a in residue2.atoms if a.name in cap_names}
    for atom in residue2.atoms:
        atom.bonds = [b for b in getattr(atom, 'bonds', [])
                      if id(b) not in remove_ids]
    residue2.atoms = [a for a in residue2.atoms if id(a) not in remove_ids]
    logger.debug("Removed cap atoms %s from fragment.", cap_names)

    # ── Step 2: rename CA candidate; assign HA ───────────────────────────────
    ca_atom = next((a for a in residue2.atoms if a.name == ca_atom_name), None)
    if ca_atom is None:
        logger.error("CA candidate '%s' not found after cap removal; "
                     "available: %s", ca_atom_name,
                     [a.name for a in residue2.atoms])
        return
    ca_atom.name = 'CA'

    ca_h_set = set(ca_H_names) if ca_H_names else set()
    ha_assigned = False
    extra_h_ids = set()
    for atom in residue2.atoms:
        if atom.name == 'CA':
            continue
        is_ca_h = (atom.name in ca_h_set) if ca_h_set else (
            atom.name.startswith('H') or
            (atom.name[0].isdigit() and len(atom.name) > 1
             and atom.name[1].upper() == 'H'))
        if not is_ca_h:
            continue
        if not ha_assigned:
            atom.name = 'HA'
            ha_assigned = True
        else:
            extra_h_ids.add(id(atom))

    if extra_h_ids:
        for atom in residue2.atoms:
            atom.bonds = [b for b in getattr(atom, 'bonds', [])
                          if id(b) not in extra_h_ids]
        residue2.atoms = [a for a in residue2.atoms if id(a) not in extra_h_ids]

    # ── Step 3: SVD rigid-body fit; clone backbone atoms ─────────────────────
    r1_pts, r2_pts = [], []
    r1_by = {a.name: a for a in residue1.atoms}
    r2_by = {a.name: a for a in residue2.atoms}

    cb1 = r1_by.get('CB')
    cb2 = r2_by.get('CB')
    if cb1 and cb2:
        r1_pts.append(np.array(cb1.x, dtype=float))
        r2_pts.append(np.array(cb2.x, dtype=float))
        cb1_hs = sorted(
            [b for b in getattr(cb1, 'bonds', [])
             if b.name.lstrip('0123456789').upper().startswith('H')],
            key=lambda a: a.name)
        cb2_hs = [r2_by[n] for n in ('HB1', 'HB2', 'HB3') if n in r2_by]
        for h1, h2 in zip(cb1_hs, cb2_hs):
            r1_pts.append(np.array(h1.x, dtype=float))
            r2_pts.append(np.array(h2.x, dtype=float))
    else:
        for name in ('CA', 'HA'):
            if name in r1_by and name in r2_by:
                r1_pts.append(np.array(r1_by[name].x, dtype=float))
                r2_pts.append(np.array(r2_by[name].x, dtype=float))

    if len(r1_pts) == 0:
        R_bb, t_bb = np.eye(3), np.zeros(3)
    elif len(r1_pts) == 1:
        R_bb, t_bb = np.eye(3), r2_pts[0] - r1_pts[0]
    else:
        old = np.array(r1_pts)
        new = np.array(r2_pts)
        oc, nc = old.mean(0), new.mean(0)
        H = (old - oc).T @ (new - nc)
        U, _, Vt = np.linalg.svd(H)
        R_bb = Vt.T @ U.T
        if np.linalg.det(R_bb) < 0:
            Vt[-1] *= -1
            R_bb = Vt.T @ U.T
        t_bb = nc - R_bb @ oc

    bb_to_clone = (_N_SIDE | {'C', 'OT1', 'OT2'}) if term == 'C' \
        else (_N_SIDE | _C_SIDE)
    existing = {a.name for a in residue2.atoms}

    n_side, c_side = [], []
    for atom in residue1.atoms:
        if atom.name not in bb_to_clone or atom.name in existing:
            continue
        cloned = deepcopy(atom)
        cloned.bonds   = []
        cloned.resname = residue2.resname
        cloned.x = list(R_bb @ np.array(atom.x, dtype=float) + t_bb)
        (n_side if atom.name in _N_SIDE else c_side).append(cloned)

    residue2.atoms = n_side + residue2.atoms + c_side
    logger.info("_prepare_nsaa_fragment: cloned %s from %s.",
                [a.name for a in n_side + c_side], residue1.resname)


# =============================================================================
# Parameter assignment
# =============================================================================

def _assign_itp_entries(residue, itp_file: str, rtp, backbone_resname=None):
    """Assign atomtype and charge to an NSAA residue from an ITP fragment.

    Backbone atoms (N, HN/H, CA, HA, C, O) receive canonical FF parameters
    from the RTP; sidechain atoms receive types/charges from the ITP.

    Returns a list of atom names that could not be assigned (empty = all OK).
    """
    itp_atoms = _parse_itp_atoms(itp_file)
    itp_bonds = _parse_itp_bonds(itp_file)
    cap_C, cap_Hs, _, _ = _find_methyl_cap(itp_atoms, itp_bonds)

    cap_set = set(cap_Hs)
    if cap_C:
        cap_set.add(cap_C)
    itp_params = {a['name']: {'type': a['type'], 'charge': a['charge'],
                               'mass': a['mass']}
                  for a in itp_atoms if a['name'] not in cap_set}

    if backbone_resname is None:
        backbone_resname = 'ALA' if 'ALA' in rtp else next(iter(rtp), None)
    bb_params = {}
    if backbone_resname and backbone_resname in rtp:
        for entry in rtp[backbone_resname]['atoms']:
            if entry[0] in _BACKBONE:
                bb_params[entry[0]] = {'type': entry[1], 'charge': float(entry[2])}

    missing = []
    for atom in residue.atoms:
        name = atom.name
        if name in bb_params:
            atom.atomtype = bb_params[name]['type']
            atom.q = bb_params[name]['charge']
        elif name in itp_params:
            atom.atomtype = itp_params[name]['type']
            atom.q = itp_params[name]['charge']
            if itp_params[name]['mass'] > 0:
                atom.m = itp_params[name]['mass']
        else:
            # Try CHARMM numeric-prefix alias (e.g. "1HB" ↔ "HB1")
            aliased = (name[1:] + name[0]) if name[0].isdigit() else None
            if aliased and aliased in itp_params:
                atom.atomtype = itp_params[aliased]['type']
                atom.q = itp_params[aliased]['charge']
                if itp_params[aliased]['mass'] > 0:
                    atom.m = itp_params[aliased]['mass']
            else:
                missing.append(name)

    if missing:
        logger.warning("No parameters found for atoms in %s: %s",
                       residue.resname, missing)
    return missing


def _assign_rtp_entries(mol, rtp, term=None):
    """Assign RTP parameters (atomtype, charge, bonds) to *mol*.

    Returns a list of inter-residue bond pairs (e.g. ``['+N', 'C']``).
    """
    resname = mol.resname
    if resname not in rtp:
        raise ValueError(
            f"Residue '{resname}' not found in RTP. "
            "Check your RTP definitions.")

    entr = _rtp_terminal(rtp[resname], term) if term else rtp[resname]
    logger.debug("RTP entry for %s: %d atoms", resname, len(entr['atoms']))

    neigh = []

    for atom_entry in entr['atoms']:
        atom_name, atom_type, atom_q = atom_entry[0], atom_entry[1], atom_entry[2]
        # Apply CHARMM alias inversion
        if resname in _inverted_aliases and atom_name in _inverted_aliases[resname]:
            atom_name = _inverted_aliases[resname][atom_name]
        elif atom_name[0].isnumeric():
            # Numeric-prefix to suffix (e.g. "1HB" → "HB1")
            atom_name = atom_name[1:] + atom_name[0]

        try:
            atom = mol.fetchm([atom_name], permissive=True)[0]
        except IndexError:
            logger.debug("Atom '%s' not found in %s — skipped.", atom_name, resname)
            continue
        atom.atomtype = atom_type
        atom.q = atom_q
        atom.cgnr = atom_entry[3]

    for a1, a2 in entr['bonds']:
        within = '-' not in a1 and '+' not in a1 and '-' not in a2 and '+' not in a2
        if within:
            if resname in _inverted_aliases:
                a1 = _inverted_aliases[resname].get(a1, a1)
                a2 = _inverted_aliases[resname].get(a2, a2)
            try:
                at1, at2 = mol.fetchm([a1, a2], permissive=True)
                at1.bonds.append(at2)
                at2.bonds.append(at1)
            except (IndexError, Exception):
                logger.debug("Bond %s--%s not found in %s.", a1, a2, resname)
        else:
            neigh.append([a1, a2])

    logger.info("Matched %s with RTP.", resname)
    return neigh


def _validate_atom_parameters(residue, context: str = '') -> List[str]:
    """Warn about atoms that lack an atomtype; return their names."""
    missing = [a.name for a in residue.atoms if not getattr(a, 'atomtype', None)]
    if missing:
        logger.warning("Atoms in %s (%s) have no atomtype: %s",
                       residue.resname, context, missing)
    return missing


def _assign_mass_atp(r1, r2, ffatomtypes: str, itp_file: str = None):
    """Assign masses to all atoms in *r1* and *r2*.

    Priority: (1) FF atomtypes.atp, (2) already-set atom.m, (3) ITP
    per-atom mass, (4) element-symbol fallback.
    """
    mass: Dict[str, float] = {}
    with open(ffatomtypes) as fp:
        for line in kickOutComments(fp.readlines(), ';'):
            parts = line.split()
            if len(parts) >= 2:
                try:
                    mass[parts[0]] = float(parts[1])
                except ValueError:
                    pass

    itp_mass = {}
    if itp_file:
        for a in _parse_itp_atoms(itp_file):
            if a['mass'] > 0:
                itp_mass[a['name']] = a['mass']

    for atom in r1.atoms + r2.atoms:
        atype = getattr(atom, 'atomtype', '')
        if atype in mass:
            atom.m = mass[atype]
        elif getattr(atom, 'm', 0) > 0:
            pass  # already set (e.g. from ITP)
        elif atom.name in itp_mass:
            atom.m = itp_mass[atom.name]
        else:
            sym = (getattr(atom, 'symbol', '') or atom.name[0]).upper()
            if sym in _ELEM_MASS:
                atom.m = _ELEM_MASS[sym]
                logger.warning("No mass for type '%s' (%s); using element "
                               "fallback %s = %.3f Da.", atype, atom.name,
                               sym, atom.m)
            else:
                logger.error("Cannot determine mass for '%s' (type='%s').",
                             atom.name, atype)


# =============================================================================
# RTP entry manipulation helpers
# =============================================================================

def drop_refs(rows: List, names) -> List:
    """Drop rows from a bond/dihedral/improper list that reference *names*."""
    names = set(names)
    return [r for r in rows if not any(x in names for x in r)]


def del_atoms(entr: dict, *names) -> dict:
    """Remove atoms (and all references to them) from an RTP entry dict."""
    names_set = set(names)
    entr['atoms'] = [a for a in entr['atoms'] if a[0] not in names_set]
    for k in ('bonds', 'diheds', 'improps', 'cmap'):
        entr[k] = drop_refs(entr.get(k, []), names_set)
    return entr


def add_atom(entr: dict, name: str, atype: str, q: float,
             after: str = None) -> dict:
    """Append or insert an atom into an RTP entry dict."""
    atoms = entr['atoms']
    idx = max(a[3] for a in atoms) + 1
    row = [name, atype, float(q), idx]
    if after is None:
        atoms.append(row)
    else:
        for i, a in enumerate(atoms):
            if a[0] == after:
                atoms.insert(i + 1, row)
                break
        else:
            raise KeyError(f"atom not found: {after}")
    return entr


def ren_atom(entr: dict, old: str, new: str = None,
             atype: str = None, q: float = None) -> dict:
    """Rename / retype / recharge an atom in an RTP entry dict."""
    for a in entr['atoms']:
        if a[0] == old:
            if new   is not None: a[0] = new
            if atype is not None: a[1] = atype
            if q     is not None: a[2] = float(q)
            break
    else:
        raise KeyError(f"atom not found: {old}")
    if new is not None and new != old:
        for k in ('bonds', 'diheds', 'improps', 'cmap'):
            for r in entr.get(k, []):
                for i, x in enumerate(r):
                    if x == old:
                        r[i] = new
    return entr


def _rtp_terminal(entr: dict, term: str) -> dict:
    """Apply CHARMM-style terminal patches to a copy of an RTP entry dict.

    Parameters
    ----------
    entr : dict
        RTP entry dict with keys ``atoms``, ``bonds``, ``diheds``,
        ``improps``, ``cmap``.
    term : {'N', 'C'}
        Which terminus to patch.

    .. note::
       This function applies **CHARMM36** terminal patches only.
       AMBER/OPLS terminal modifications are different and not supported here.
    """
    entr = deepcopy(entr)
    if term == 'N':
        ren_atom(entr, 'N',  atype='NH3', q=-0.30)
        ren_atom(entr, 'CA', atype='CT1', q=0.21)
        ren_atom(entr, 'HA', atype='HB1', q=0.10)
        del_atoms(entr, 'HN')
        del_atoms(entr, '-C')
        add_atom(entr, 'H1', 'HC', 0.33, after='N')
        add_atom(entr, 'H2', 'HC', 0.33, after='H1')
        add_atom(entr, 'H3', 'HC', 0.33, after='H2')
        entr['bonds'].extend([['N', 'H1'], ['N', 'H2'], ['N', 'H3']])
    elif term == 'C':
        ren_atom(entr, 'C', atype='CC',  q=0.34)
        ren_atom(entr, 'O', new='OT1', atype='OC', q=-0.67)
        add_atom(entr, 'OT2', 'OC', -0.67, after='OT1')
        del_atoms(entr, '+N')
        entr['bonds'].extend([['C', 'OT1'], ['C', 'OT2']])
        entr['improps'].append(['C', 'CA', 'OT2', 'OT1', ''])
    return entr


# =============================================================================
# Atom mapping (requires RDKit)
# =============================================================================

def _reformat_pdb(mdl, strict: bool = False) -> Tuple[str, dict, list]:
    """Return a PDB-format string, original atom name map, and sigma-hole IDs."""
    m = deepcopy(mdl)
    atom_name_id: Dict[int, str] = {}
    sigma_hole_ids = []
    sh_counter = 1
    for a in m.atoms:
        if 'EP' in a.name or 'LPH' in a.name:
            a.name = f'HSH{sh_counter}'
            sh_counter += 1
            sigma_hole_ids.append(a.id)
        atom_name_id[a.id] = a.name

    at_len = 3 if strict else 4
    lines = []
    for atom in m.atoms:
        name = atom.name.upper()
        if 'CL' in name:
            atom.name = 'CL'
        elif 'BR' in name:
            atom.name = 'BR'
        elif len(atom.name) > at_len:
            atom.name = atom.name[:at_len]
        lines.append(str(atom))
    lines.append('TER\nENDMDL\n')
    return '\r\n'.join(lines), atom_name_id, sigma_hole_ids


def _set_backbone_carbonyl_double(mol, ca_idx: int, c_idx: int) -> bool:
    """Force the C=O bond type to DOUBLE in an RDKit mol."""
    c = mol.GetAtomWithIdx(c_idx)
    o_nbrs = [n for n in c.GetNeighbors() if n.GetSymbol() == 'O']
    if len(o_nbrs) >= 2:
        return False
    if len(o_nbrs) == 1:
        b = mol.GetBondBetweenAtoms(c_idx, o_nbrs[0].GetIdx())
        if b:
            b.SetBondType(Chem.BondType.DOUBLE)
            return True
    return False


def _force_backbone_pairs(n1: list, n2: list,
                          backbone1: list, backbone2: list):
    """Enforce that backbone1[k] maps to backbone2[k] in the atom mapping."""
    bb12 = dict(zip(backbone1, backbone2))
    bb21 = dict(zip(backbone2, backbone1))

    m12: Dict[int, int] = {}
    used2: set = set()
    for i1, i2 in zip(n1, n2):
        if i1 in bb12 and bb12[i1] != i2:
            continue
        if i2 in bb21 and bb21[i2] != i1:
            continue
        if i2 in used2:
            continue
        m12[i1] = i2
        used2.add(i2)

    for i1, i2 in zip(backbone1, backbone2):
        for k, v in list(m12.items()):
            if v == i2 and k != i1:
                del m12[k]
        m12[i1] = i2

    n1_out = list(m12.keys())
    return n1_out, [m12[i] for i in n1_out]


def domapping(res1, res2, backbone1: list, backbone2: list, term=None):
    """Compute atom mapping between two residues using MCS + 3-D alignment.

    Uses RDKit for MCS, alignment, and conformer embedding.

    Parameters
    ----------
    res1, res2 : Residue-like objects with ``.atoms`` list
    backbone1, backbone2 : lists of int
        Paired atom indices that must map to each other (backbone constraints).
    term : {'N', 'C', None}

    Returns
    -------
    n1, n2 : lists of int
        Paired atom indices in res1 and res2.
    pdb1, pdb2 : str
        PDB block strings for the aligned residues.
    """
    if not _RDKIT_AVAILABLE:
        raise RuntimeError(
            "domapping requires RDKit. Install it with: "
            "conda install -c conda-forge rdkit")

    res1_block, _, sigma1 = _reformat_pdb(res1)
    res2_block, _, sigma2 = _reformat_pdb(res2)

    mol1 = Chem.MolFromPDBBlock(res1_block, removeHs=False, sanitize=True)
    mol2 = Chem.MolFromPDBBlock(res2_block, removeHs=False, sanitize=True)

    _set_backbone_carbonyl_double(mol1, backbone1[1], backbone1[2])
    _set_backbone_carbonyl_double(mol2, backbone2[1], backbone2[2])

    for mol in (mol1, mol2):
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            pass
        try:
            rdmolops.AssignAtomChiralTagsFromStructure(mol)
            rdmolops.AssignStereochemistry(mol)
        except Exception:
            pass

    og_mol1 = deepcopy(mol1)
    og_mol2 = deepcopy(mol2)

    has_rings = any(a.IsInRing() for a in mol1.GetAtoms()) and \
                any(a.IsInRing() for a in mol2.GetAtoms())

    with open(os.devnull, 'w') as dev_null:
        lig_map = LigandAtomMapping(
            mol1=mol1, mol2=mol2,
            molForMcs1=deepcopy(mol1), molForMcs2=deepcopy(mol2),
            bH2H=True, bH2Hpolar=False, bH2Heavy=False,
            bdMCS=False, bRingsOnly=False, d=0.05, bChiral=True,
            sigmaHoleID1=sigma1, sigmaHoleID2=sigma2,
            timeout=10, logfile=dev_null,
            bElements=False, bCarbonize=False,
            commandName='atomMapping',
            backbone1=backbone1, backbone2=backbone2)

        # MCS pass A: element-matched, modified mols
        lig_map.bElements, lig_map.bCarbonize, lig_map.bRingsOnly = True, False, False
        n1A, n2A = lig_map.mcs()

        # MCS pass B: element-matched, original mols
        lig_map.bCarbonize = False
        n1B, n2B = lig_map.mcs()

        n1mcs, n2mcs = lig_map._compare_mappings_by_size(mol1, mol2, n1A, n2A, n1B, n2B)

        if has_rings:
            lig_map.bElements, lig_map.bCarbonize, lig_map.bRingsOnly = False, False, True
            n1C, n2C = lig_map.mcs()
            n1mcs, n2mcs = lig_map._compare_mappings_by_size(
                mol1, mol2, n1mcs, n2mcs, n1C, n2C)

        bMCSfailed = (len(n1mcs) == 0)

        n1align, n2align = [], []
        if bMCSfailed:
            lig_map.bRingsOnly = False
            n1align, n2align, _ = lig_map.alignment()

    n1mcs,  n2mcs  = _force_backbone_pairs(n1mcs,  n2mcs,  backbone1, backbone2)
    n1align, n2align = _force_backbone_pairs(n1align, n2align, backbone1, backbone2)

    if len(n1align) >= len(n1mcs):
        n1, n2 = n1align, n2align
        logger.info("Alignment mapping selected (size %d).", len(n1))
    else:
        n1, n2 = n1mcs, n2mcs
        logger.info("MCS mapping selected (size %d).", len(n1))

    if len(n1) != len(n2):
        raise RuntimeError("Mapping mismatch: |n1|=%d != |n2|=%d" % (len(n1), len(n2)))

    # Embed 100 conformers and pick the one that best aligns to res1
    cids = list(EmbedMultipleConfs(og_mol2, numConfs=100))
    if not cids:
        raise ValueError("Could not embed conformers for the target residue.")

    weights = [3 / len(n2) if (j in backbone2 and i in backbone1) else 1 / len(n2)
               for i, j in zip(n1, n2)]
    best_cid, best_rms = 0, 1e10
    for cid in cids:
        rms = AlignMol(og_mol2, og_mol1, prbCid=cid,
                       atomMap=list(zip(n2, n1)), weights=weights)
        if rms < best_rms:
            best_rms, best_cid = rms, cid

    og_mol2_best = Chem.Mol(og_mol2, confId=best_cid)
    return n1, n2, Chem.MolToPDBBlock(og_mol1), Chem.MolToPDBBlock(og_mol2_best)


# =============================================================================
# Topology / hybrid construction helpers
# =============================================================================

def _rename_model_charmm(m) -> None:
    """Rename residues and atoms to CHARMM conventions in-place."""
    rename = {'HIE': 'HSE', 'HID': 'HSD', 'HIP': 'HSP',
              'ASH': 'ASPP', 'GLH': 'GLUP', 'LYN': 'LSN'}
    for res in m.residues:
        if res.resname in rename:
            res.resname = rename[res.resname]
            for a in res.atoms:
                a.resname = res.resname
    for atom in m.atoms:
        if atom.name == 'H':
            atom.name = 'HN'
        if atom.name == 'HG' and atom.resname in ('CYS', 'SER'):
            atom.name = 'HG1'


def _merge_molecules(r1, dummies: list) -> None:
    """Append dummy atoms (B-state-only atoms from *dummies*) onto *r1*."""
    dum_h_counter = 1
    for atom in dummies:
        new = atom.copy()
        new.atomtypeB = new.atomtype
        new.qB        = new.q
        new.mB        = new.m
        new.typeB     = new.type
        new.atomtype  = 'DUM_' + new.atomtype
        new.q         = 0.0
        new.nameB     = new.name
        if new.name.startswith('H') or (new.name[0].isdigit() and
                                         len(new.name) > 1 and
                                         new.name[1] == 'H'):
            new.name = 'HV' + str(dum_h_counter)
            dum_h_counter += 1
        elif len(new.name) == 4:
            new.name = 'D' + new.name[:3]
            if new.name[1].isdigit():
                new.name = new.name[0] + new.name[2:] + new.name[1]
        else:
            new.name = 'D' + atom.name
        r1.append(new)
        logger.debug("dummy: %s -> %s", atom.name, new.name)


def _make_bstate_dummies(r1) -> None:
    """Set B-state to dummy for atoms that have not been assigned a B-state."""
    for atom in r1.atoms:
        if not hasattr(atom, 'nameB'):
            atom.nameB    = atom.name + '.gone'
            atom.atomtypeB = 'DUM_' + atom.atomtype
            atom.qB       = 0.0
            atom.mB       = atom.m


def _make_transition_dics(atom_pairs: list, r1):
    """Build A→B and B→A name dictionaries from mapped atom pairs."""
    abdic, badic = {}, {}
    for a1, a2 in atom_pairs:
        abdic[a1.name] = a2.name
        badic[a2.name] = a1.name
    for atom in r1.atoms:
        if atom.name.startswith('D') or atom.name.startswith('HV'):
            abdic[atom.name] = atom.name
            badic[atom.name] = atom.name
    return abdic, badic


def _find_atom_by_nameB(r, name: str):
    """Return the first atom in *r* whose ``nameB`` equals *name*, or None."""
    return next((a for a in r.atoms if a.nameB == name), None)


def _update_bond_lists(r1, badic: dict) -> None:
    """Reroute bond lists of dummy atoms through the B→A name map."""
    for atom in r1.atoms:
        if not (atom.name.startswith('D') or atom.name.startswith('HV')):
            continue
        new_list = []
        while atom.bonds:
            at = atom.bonds.pop(0)
            logger.debug('%s -> %s', atom.name, at.name)
            if at.name in badic:
                aa = r1.fetch(badic[at.name])[0]
                new_list.append(aa)
            else:
                aa = _find_atom_by_nameB(r1, at.name)
                if aa is not None:
                    new_list.append(aa)
                else:
                    raise RuntimeError(
                        f"Atom not found for bond rewiring: {at.name}")
        atom.bonds = new_list
        for at in atom.bonds:
            if atom not in at.bonds:
                at.bonds.append(atom)


# =============================================================================
# Bonded-term entry generation
# =============================================================================

def _improps_as_atoms(im: list, r, r2, use_b: bool = False) -> list:
    """Convert improper/dihedral name lists to lists of Atom objects."""
    result = []
    for ii in im:
        atoms_in = ii[:4]
        new_ii = []
        for name in atoms_in:
            if name and name[0] in ('+', '-'):
                new_ii.append(Atom(name=name))
                continue
            if use_b:
                if r2.resname in _inverted_aliases:
                    name = _inverted_aliases[r2.resname].get(name, name)
                a = next((at for at in r.atoms if at.nameB == name), None)
                if a is None:
                    logger.debug("Improper B-atom '%s' not found.", name)
                    a = Atom(name=name)
            else:
                if r.resname in _inverted_aliases:
                    name = _inverted_aliases[r.resname].get(name, name)
                found = r.fetch(name)
                a = found[0] if found else Atom(name=name)
            new_ii.append(a)
        new_ii.extend(ii[4:])
        result.append(new_ii)
    return result


def _improp_entries_match(lst1: list, lst2: list) -> bool:
    return (all(a.name == b.name for a, b in zip(lst1, lst2)) or
            all(a.name == b.name for a, b in zip(lst1, reversed(lst2))))


def _generate_dihedral_entries(im1: list, im2: list, r, pairs: list) -> list:
    """Merge A-state and B-state dihedral lists into hybrid entries."""
    result, done1, done2 = [], [], []
    for i1 in im1:
        for i2 in im2:
            if _improp_entries_match(i1[:4], i2[:4]) and i2 not in done2:
                result.append(i1[:4] + [i1[4] or 'default-A', i2[4] or 'default-B'])
                done1.append(i1)
                done2.append(i2)
                break
    for i1 in im1:
        if i1 in done1:
            continue
        entry = i1[:4] + [i1[4] or 'default-A']
        gone = any('gone' in getattr(a, 'nameB', '') for a in i1[:4])
        if gone:
            entry.append(i1[4] or 'default-A')
        elif 'torsion' in (i1[4] or ''):
            entry.append('un' + i1[4].replace('torsion', 'tors'))
        elif 'dih_' in (i1[4] or ''):
            entry.append('un' + i1[4])
        else:
            entry.append('un')
        result.append(entry)
    for i2 in im2:
        if i2 in done2:
            continue
        special = any(i2[x].name.startswith(('D', 'HV')) for x in range(4))
        if i2[4] == '':
            a_label = 'default-B' if special else 'un'
            result.append(i2[:4] + [a_label, 'default-B'])
        else:
            b_label = i2[4] if special else 'un' if 'torsion' not in i2[4] and 'dih_' not in i2[4] else 'un' + i2[4]
            result.append(i2[:4] + [i2[4], b_label])
    return result


def _generate_improp_entries(im1: list, im2: list, r) -> list:
    """Merge A-state and B-state improper lists into hybrid entries."""
    result, done1, done2 = [], [], []
    for i1 in im1:
        for i2 in im2:
            if _improp_entries_match(i1[:4], i2[:4]) and i2 not in done2:
                a_lbl = 'default-star' if i1[4] == '105.4' else (i1[4] or 'default-A')
                b_lbl = 'default-star' if i2[4] == '105.4' else (i2[4] or 'default-B')
                result.append(i1[:4] + [a_lbl, b_lbl])
                done1.append(i1)
                done2.append(i2)
                break
    for i1 in im1:
        if i1 in done1:
            continue
        gone = any(hasattr(a, 'nameB') and 'gone' in a.nameB for a in i1[:4])
        if i1[4] == '':
            result.append(i1[:4] + ['default-A', 'default-A' if gone else 'un'])
        elif i1[4] == '105.4':
            result.append(i1[:4] + ['default-star', 'un'])
        else:
            result.append(i1[:4] + [i1[4], i1[4] if gone else 'un'])
    for i2 in im2:
        if i2 in done2:
            continue
        special = any(i2[x].name.startswith(('D', 'HV')) for x in range(4))
        if i2[4] == '':
            result.append(i2[:4] + ['default-B' if special else 'un', 'default-B'])
        elif i2[4] == '105.4':
            result.append(i2[:4] + ['un', 'default-star'])
        else:
            result.append(i2[:4] + [i2[4] if special else 'un', i2[4]])
    return result


# =============================================================================
# File generation
# =============================================================================

def _make_rtp(r, ii_list: list, dihi_list: list,
              neigh_bonds: list, cmap: list, ff: str) -> List[str]:
    """Generate RTP file content for the hybrid residue *r*."""
    header = _CHARMM_RTP_HEADER if 'charmm' in ff else _AMBER_RTP_HEADER
    rtp = [header, f'\n[ {r.resname} ] ; {r.resnA} -> {r.resnB}\n\n']

    rtp.append('[ atoms ]\n')
    for cgnr, atom in enumerate(r.atoms, 1):
        rtp.append(f'   {atom.name:6}   {atom.atomtype:<15}  {atom.q:8.5f}  {cgnr}\n')

    rtp.append('\n[ bonds ]\n')
    for atom in r.atoms:
        for at in atom.bonds:
            if atom.id < at.id:
                rtp.append(f'   {atom.name:6}  {at.name:6}'
                           f' ; ({atom.nameB:10}  {at.nameB:<10})\n')
    for pair in neigh_bonds:
        rtp.append(f'   {pair[0]:6}  {pair[1]:6}\n')

    def _fmt(section, interactions):
        rtp.append(f'\n[ {section} ]\n')
        for ii in interactions:
            if not ii[4].startswith('default'):
                rtp.append(f'   {ii[0].name:6}  {ii[1].name:6}  '
                           f'{ii[2].name:6}  {ii[3].name:6}  {ii[4]}\n')
            else:
                rtp.append(f'   {ii[0].name:6}  {ii[1].name:6}  '
                           f'{ii[2].name:6}  {ii[3].name:6}\n')

    _fmt('impropers', ii_list)
    _fmt('dihedrals', dihi_list)

    if cmap:
        rtp.append('\n[ cmap ]\n')
        for item in cmap:
            rtp.append(f'   {item}\n')
    return rtp


def _make_mtp(r, ii_list: list, dihi_list: list,
              rotations=None) -> List[str]:
    """Generate MTP file content for the hybrid residue *r*."""
    mtp = [f'\n[ {r.resname} ] ; {r.resnA} -> {r.resnB}\n\n']

    mtp.append('\n[ morphes ]\n')
    for atom in r.atoms:
        mtp.append(f'   {atom.name:10} {atom.atomtype:10} -> '
                   f'{atom.nameB:10} {atom.atomtypeB:10}\n')

    mtp.append('\n[ atoms ]\n')
    for cgnr, atom in enumerate(r.atoms, 1):
        tag = (' ; types != | charge !='
               if atom.atomtype != atom.atomtypeB or atom.q != atom.qB
               else ' ; types == | charge ==')
        mtp.append(f'   {atom.name:8} {atom.atomtype:10} {atom.q:10.6f} '
                   f'{cgnr:6} {atom.m:10.6f} {atom.atomtypeB:10} '
                   f'{atom.qB:10.6f} {atom.mB:10.6f}  {tag}\n')

    mtp.append('\n[ coords ]\n')
    for atom in r.atoms:
        mtp.append(f'   {atom.x[0]:8.3f} {atom.x[1]:8.3f} {atom.x[2]:8.3f}\n')

    def _fmt(section, interactions):
        mtp.append(f'\n[ {section} ]\n')
        for ii in interactions:
            mtp.append(f' {ii[0].name:6} {ii[1].name:6} {ii[2].name:6} '
                       f'{ii[3].name:6}     {ii[4]:25} {ii[5]:25}\n')

    _fmt('impropers', ii_list)
    _fmt('dihedrals', dihi_list)

    if rotations:
        mtp.append('\n[ rotations ]\n')
        for rot in rotations:
            mtp.append(f'  {rot[0].name}-{rot[1].name} '
                       f'{" ".join(a.name for a in rot[2:])}\n')
    return mtp


def _get_residue_types(filename: str) -> set:
    """Return the set of residue names listed in *filename* (e.g. residuetypes.dat)."""
    types = set()
    if os.path.isfile(filename):
        with open(filename) as f:
            for line in f:
                line = line.lstrip()
                if line and not line.startswith((';', '#')):
                    types.add(line.split()[0])
    return types


def _get_existing_atomtypes(filename: str) -> set:
    """Return atom types present in the ``[ atomtypes ]`` section of *filename*."""
    types = set()
    in_section = False
    if os.path.isfile(filename):
        with open(filename) as f:
            for line in f:
                s = line.strip()
                if s.startswith('['):
                    in_section = (s == '[ atomtypes ]')
                elif in_section and s and not s.startswith((';', '#')):
                    parts = s.split()
                    if parts:
                        types.add(parts[0])
    return types


def _write_atp_fnb(fn_atp: str, fn_nb: str, r, ffpath: str) -> None:
    """Append new DUM atom types to the FF ``atomtypes.atp`` and
    ``ffnonbonded.itp`` files.  Existing types are not duplicated.
    """
    atp_path = os.path.join(ffpath, fn_atp)
    nb_path  = os.path.join(ffpath, fn_nb)

    existing_atp = _get_residue_types(atp_path)
    existing_nb  = _get_existing_atomtypes(nb_path)

    # atomtypes.atp
    with open(atp_path, 'a' if os.path.isfile(atp_path) else 'w') as f:
        for atom in r.atoms:
            for atype, mass in ((atom.atomtype, atom.m),
                                (atom.atomtypeB, atom.mB)):
                if atype.startswith('DUM') and atype not in existing_atp:
                    f.write(f'{atype:<6}  {mass:10.6f}\n')
                    existing_atp.add(atype)

    # ffnonbonded.itp — insert into [ atomtypes ] section
    existing_lines = []
    if os.path.isfile(nb_path):
        with open(nb_path) as f:
            existing_lines = f.readlines()

    modified = list(existing_lines)
    in_at = False
    insert_idx = None
    for i, line in enumerate(modified):
        s = line.strip()
        if s.startswith('['):
            in_at = (s == '[ atomtypes ]')
            if in_at:
                insert_idx = i + 1
        elif in_at and s and not s.startswith((';', '#')):
            parts = s.split()
            if parts:
                existing_nb.add(parts[0])

    new_lines = []
    if insert_idx is not None:
        for atom in r.atoms:
            for atype, mass in ((atom.atomtype, atom.m),
                                (atom.atomtypeB, atom.mB)):
                if atype.startswith('DUM') and atype not in existing_nb \
                        and atype not in [l.split()[0] for l in new_lines]:
                    new_lines.append(
                        f'{atype:<10}\t0\t{mass:4.2f}\t   0.0000  A'
                        f'   0.00000e+00 0.00000e+00\n')
                    existing_nb.add(atype)
        for j, line in enumerate(new_lines):
            modified.insert(insert_idx + j, line)

    with open(nb_path, 'w') as f:
        f.writelines(modified)


# =============================================================================
# Main entry point
# =============================================================================

def mutate_nsaa(m, mut_resid: int, mut_resname: str, ff: str,
                itp: str = None, nonstandardPDB: str = None,
                mut_chain: str = None, inplace: bool = False,
                verbose: bool = False) -> None:
    """Mutate residue *mut_resid* to a non-standard amino acid.

    Generates a hybrid structure and writes the following files:

    * ``<hybrid>.rtp`` — into the FF directory (read by pdb2gmx)
    * ``<hybrid>.mtp`` — in the working directory (pass to gen_hybrid_top)
    * ``<hybrid>.arn`` — atom-renaming spec for CHARMM FFs (into FF dir)

    Updates ``atomtypes.atp`` and ``ffnonbonded.itp`` in the FF directory
    with any new DUM atom types.

    Parameters
    ----------
    m : Model
        Input protein model.  Must have been loaded with ``rename_atoms=True``
        and ``scale_coords='A'``.
    mut_resid : int
        Residue ID to mutate (after renumbering, if applicable).
    mut_resname : str
        Name of the target NSAA (e.g. ``'SEP'``, ``'TPO'``).
    ff : str
        Force-field name (e.g. ``'charmm36m-mut'``).
    itp : str, optional
        Path to a GROMACS ITP fragment file for the NSAA, parameterized with
        a CH3 cap at the CA position.
    nonstandardPDB : str, optional
        Path to a PDB/GRO file for the NSAA.  If omitted, the NSAA is built
        from the pmx library (standard residues only).
    mut_chain : str, optional
        Chain ID; required only if the model has multiple chains.
    inplace : bool
        Modify *m* in place.  Default is False (returns nothing either way
        since the mutation is a side effect on the structure).
    verbose : bool
        Extra logging.
    """
    if not _RDKIT_AVAILABLE:
        raise RuntimeError(
            "mutate_nsaa requires RDKit. Install with:\n"
            "    conda install -c conda-forge rdkit")

    from .builder import build_chain

    logger.info("Starting NSAA mutation: res %d -> %s", mut_resid, mut_resname)

    ffpath = get_ff_path(ff)
    bCharmm = 'charmm' in ff

    m2 = m if inplace else deepcopy(m)
    if m2.unity == 'nm':
        m2.nm2a()

    residue1 = m2.fetch_residue(idx=mut_resid, chain=mut_chain)
    chain    = residue1.chain
    cterm    = chain.cterminus()
    nterm    = chain.nterminus()
    term     = 'C' if residue1 == cterm else ('N' if residue1 == nterm else None)

    # ── Load the target residue ───────────────────────────────────────────────
    model2 = Model()
    if nonstandardPDB is None:
        if mut_resname.upper() not in library._aacids_dic.values():
            raise ValueError(
                f"No structure file provided for NSAA '{mut_resname}' and it "
                "is not in the pmx standard library.")
        one = library._one_letter.get(mut_resname, mut_resname)
        ch  = build_chain(one)
        if term == 'C':
            ch.add_cterm_charged()
        model2.residues = ch.residues
        model2.moltype  = 'protein'
        model2.atoms    = model2.residues[0].atoms
        if bCharmm:
            _rename_model_charmm(model2)
    else:
        model2.read(nonstandardPDB)
    residue2 = model2.residues[0]

    # ── Parse force-field RTP ─────────────────────────────────────────────────
    try:
        rtp = __import__('pmx.ffparser', fromlist=['RTPParser']).RTPParser(
            os.path.join(ffpath, 'merged.rtp'))
    except Exception:
        from .ffparser import RTPParser
        rtp = RTPParser(os.path.join(ffpath, 'aminoacids.rtp'))

    bond_neigh = _assign_rtp_entries(residue1, rtp, term)

    # ── Assign parameters to the target residue ───────────────────────────────
    if residue2.resname in rtp:
        logger.info("'%s' found in RTP — using RTP parameters.", residue2.resname)
        try:
            _assign_rtp_entries(residue2, rtp, term)
        except Exception as e:
            logger.warning("RTP assignment for %s failed: %s", residue2.resname, e)
    elif itp is not None:
        logger.info("'%s' not in RTP — using ITP fragment: %s",
                    residue2.resname, itp)
        itp_atoms_raw = _parse_itp_atoms(itp)
        itp_bonds_raw = _parse_itp_bonds(itp)
        cap_C, cap_Hs, cap_CA, cap_CA_Hs = _find_methyl_cap(itp_atoms_raw,
                                                             itp_bonds_raw)
        has_bb = all(any(a.name == n for a in residue2.atoms)
                     for n in ('N', 'C', 'O'))
        if not has_bb:
            if cap_C is None:
                raise RuntimeError(
                    "Capped-fragment workflow: residue2 has no backbone N/C/O "
                    "and no CH3 cap was detected in the ITP.")
            _prepare_nsaa_fragment(residue2, residue1,
                                   cap_C, cap_Hs,
                                   ca_atom_name=cap_CA,
                                   ca_H_names=cap_CA_Hs,
                                   term=term)
        _assign_itp_entries(residue2, itp, rtp,
                            backbone_resname=residue1.resname)
    else:
        logger.warning("'%s' not in RTP and no ITP provided — "
                       "parameters may be incomplete.", residue2.resname)

    _validate_atom_parameters(residue1, 'state A')
    _validate_atom_parameters(residue2, 'state B')
    _assign_mass_atp(residue1, residue2,
                     os.path.join(ffpath, 'atomtypes.atp'), itp_file=itp)

    # ── Backbone index lists for atom mapping ─────────────────────────────────
    bb_names = ['N', 'CA', 'C']
    bb_names += (['OT1', 'OT2'] if term == 'C' else ['O'])
    if (any(a.name == 'CB' for a in residue1.atoms) and
            any(a.name == 'CB' for a in residue2.atoms)):
        bb_names.append('CB')

    backbone1, backbone2 = [], []
    for name in bb_names:
        i1 = next((c for c, a in enumerate(residue1.atoms) if a.name == name), None)
        i2 = next((c for c, a in enumerate(residue2.atoms) if a.name == name), None)
        if i1 is not None and i2 is not None:
            backbone1.append(i1)
            backbone2.append(i2)
    logger.info("Backbone mapping: %s", list(zip(backbone1, backbone2)))

    # ── MCS atom mapping ──────────────────────────────────────────────────────
    pairs1, pairs2, res1_aligned, res2_aligned = \
        domapping(residue1, residue2, backbone1, backbone2, term)

    res3 = Model(pdbline=res1_aligned).residues[0]
    res4 = Model(pdbline=res2_aligned).residues[0]
    for c, at in enumerate(residue1.atoms):
        at.x = res3.atoms[c].x
    for c, at in enumerate(residue2.atoms):
        at.x = res4.atoms[c].x

    # Ensure polar H is in the mapping
    hn_name = 'HN' if bCharmm else 'H'
    for c, a in enumerate(residue1.atoms):
        if a.name == hn_name:
            pairs1.append(c)
            break
    for c, a in enumerate(residue2.atoms):
        if a.name == hn_name:
            pairs2.append(c)
            break

    # ── Build hybrid residue ──────────────────────────────────────────────────
    merged_atoms2 = []
    residue1.batoms = []
    atom_pairs = []

    for id1, id2 in zip(pairs1, pairs2):
        at1, at2 = residue1.atoms[id1], residue2.atoms[id2]
        logger.info("  %s ---> %s", at1.name, at2.name)
        at1.atomtypeB, at1.qB, at1.mB, at1.nameB = \
            at2.atomtype, at2.q, at2.m, at2.name
        residue1.batoms.append(at2)
        merged_atoms2.append(at2)
        atom_pairs.append([at1, at2])

    dummies = [a for a in residue2.atoms if a not in merged_atoms2]
    logger.info("Dummy atoms (B-state only): %s", [a.name for a in dummies])
    _merge_molecules(residue1, dummies)
    _make_bstate_dummies(residue1)

    _write_atp_fnb('atomtypes.atp', 'ffnonbonded.itp', residue1, ffpath)

    # ── Generate bonded-term entries ──────────────────────────────────────────
    abdic, badic = _make_transition_dics(atom_pairs, residue1)
    _update_bond_lists(residue1, badic)

    # Use _rtp_terminal for NSAA target only if it's in RTP
    entr1 = _rtp_terminal(rtp[residue1.resname], term) if term \
        else rtp[residue1.resname]
    if residue2.resname in rtp:
        entr2 = _rtp_terminal(rtp[residue2.resname], term) if term \
            else rtp[residue2.resname]
    else:
        entr2 = {'atoms': [], 'bonds': [], 'diheds': [], 'improps': [], 'cmap': []}

    dih1 = _improps_as_atoms(entr1['diheds'], residue1, residue2)
    dih2 = _improps_as_atoms(entr2['diheds'], residue1, residue2, use_b=True)
    im1  = _improps_as_atoms(entr1['improps'], residue1, residue2)
    im2  = _improps_as_atoms(entr2['improps'], residue1, residue2, use_b=True)

    cmap      = rtp[residue1.resname].get('cmap', []) if bCharmm else []
    dihi_list = _generate_dihedral_entries(dih1, dih2, residue1, atom_pairs)
    ii_list   = _generate_improp_entries(im1, im2, residue1)

    # ── Name and write outputs ────────────────────────────────────────────────
    one_letter = library._one_letter.get(residue1.resname, residue1.resname[0])
    hybrid_name = one_letter + '2' + mut_resname[:3].upper()

    residue1.resnA, residue1.resnB = residue1.resname, mut_resname
    residue1.set_resname(hybrid_name)

    rtp_content = _make_rtp(residue1, ii_list, dihi_list, bond_neigh, cmap, ff)
    mtp_content = _make_mtp(residue1, ii_list, dihi_list, rotations=None)

    list2file(rtp_content, os.path.join(ffpath, f'{hybrid_name}.rtp'))
    list2file(mtp_content, f'./{hybrid_name}.mtp')
    if bCharmm:
        list2file([_ARN_CONTENT], os.path.join(ffpath, f'{hybrid_name}.arn'))

    # Register the hybrid residue name in residuetypes.dat
    fn_restypes = pmx_data_file('mutff/residuetypes.dat')
    existing    = _get_residue_types(fn_restypes)
    mode = 'a' if os.path.isfile(fn_restypes) else 'w'
    with open(fn_restypes, mode) as f:
        if hybrid_name not in existing:
            f.write(f'{hybrid_name:<6} Protein\n')

    logger.info("NSAA mutation complete: hybrid residue '%s' written.", hybrid_name)
    logger.info("  RTP  -> %s/%s.rtp", ffpath, hybrid_name)
    logger.info("  MTP  -> ./%s.mtp  (pass to gen_hybrid_top)", hybrid_name)


# =============================================================================
# ITP → FF conversion  (RDKit-free path)
# =============================================================================

# AMBER → GAFF type mapping for the most common sidechain atom types.
# Used by sanitize_itp to remap abstract OpenFF types for atoms that
# correspond to a known sidechain position in the base residue.
_AMBER_TO_GAFF: Dict[str, str] = {
    # sp3 carbons
    'CT': 'c3', 'CT1': 'c3', 'CT2': 'c3', 'CT3': 'c3', 'CX': 'c3',
    # sp2 carbons
    'C':  'c',  'C*': 'cc', 'C5': 'cc', 'C4': 'cd',
    # sp3 nitrogens
    'N3': 'n3', 'NT': 'n3',
    # amide nitrogens
    'NH1': 'n', 'NH2': 'n',
    # hydroxyl oxygen
    'OH1': 'oh', 'OH': 'oh',
    # thiol sulfur
    'S':  'sh', 'SH1': 'sh',
    # thioether sulfur
    'S1': 'ss',
    # polar H
    'H':  'h1', 'HO': 'ho', 'HS': 'hs',
    # aliphatic H
    'HA': 'hc', 'HB': 'hc', 'HP': 'hc',
}


def sanitize_itp(
    itp_file: str,
    base_resname: str,
    ff: str,
    atomtypes_itp: Optional[str] = None,
    cap_C_override: Optional[str] = None,
    out_itp: Optional[str] = None,
) -> Tuple[str, dict]:
    """Normalize atom names and types in an OpenFF or CGenFF ITP fragment.

    Two problems are fixed:

    1. **Non-unique atom names** — OpenFF/CGenFF often uses generic element
       symbols (``C``, ``N``, ``H``) for every atom.  These are renamed to
       element+nr (``C5``, ``N6``, ``H23``) so that cap detection and bonded
       wiring are unambiguous.

    2. **Abstract atom types** — OpenFF uses ``AT_0``, ``AT_1`` … that carry
       no AMBER/GAFF meaning.  Sidechain atoms that correspond to a known
       base-residue position get the GAFF equivalent via :data:`_AMBER_TO_GAFF`.
       Pure-ligand atoms (beyond the base sidechain) get the closest GAFF type
       by sigma/epsilon nearest-neighbour matching (requires *atomtypes_itp*).

    Parameters
    ----------
    itp_file : str
        Path to the raw ITP fragment.
    base_resname : str
        Base amino acid (e.g. ``'CYS'``).
    ff : str
        Force-field name.
    atomtypes_itp : str or None
        Path to the OpenFF/CGenFF ``[ atomtypes ]`` ITP.  Required for
        sigma/epsilon-based type remapping of pure-ligand atoms; ignored
        for GAFF inputs.
    cap_C_override : str or None
        Override the auto-detected cap carbon name.
    out_itp : str or None
        Output path (default: ``<stem>_sanitized.itp``).

    Returns
    -------
    out_itp_path : str
    itp_map : dict
        ``{atom_nr: {'new_name': str, 'amber_type': str or None,
        'gaff_type': str}}``
    """
    from collections import deque, defaultdict as _dd
    from .ffparser import RTPParser

    if out_itp is None:
        base, ext = os.path.splitext(itp_file)
        out_itp = base + '_sanitized' + ext

    # ── 1. Parse ITP ────────────────────────────────────────────────────────
    itp_atoms    = _parse_itp_atoms(itp_file)
    _make_itp_names_unique(itp_atoms)
    itp_bonds_raw = _parse_itp_bonds(itp_file)

    nr_to_atom = {a['nr']: a for a in itp_atoms}
    adj_nr: Dict[int, set] = {a['nr']: set() for a in itp_atoms}
    for ai, aj in itp_bonds_raw:
        adj_nr[ai].add(aj)
        adj_nr[aj].add(ai)

    # ── 2. Find CH3 cap ──────────────────────────────────────────────────────
    cap_C_name, cap_H_names, connect_name, _ = _find_methyl_cap(
        itp_atoms, itp_bonds_raw)
    if cap_C_override:
        cap_C_name = cap_C_override
        cap_C_nr   = next((a['nr'] for a in itp_atoms if a['name'] == cap_C_name), None)
        cap_H_names = [nr_to_atom[n]['name'] for n in adj_nr.get(cap_C_nr, set())
                       if nr_to_atom[n]['name'].startswith('H')]
        non_h = [n for n in adj_nr.get(cap_C_nr, set())
                 if not nr_to_atom[n]['name'].startswith('H')]
        connect_name = nr_to_atom[non_h[0]]['name'] if non_h else connect_name

    if cap_C_name is None:
        raise ValueError(
            f"Cannot auto-detect CH3 cap in '{itp_file}'. "
            "Use cap_C_override to specify the cap carbon.")

    cap_all_names = {cap_C_name} | set(cap_H_names or [])
    cap_nrs = {a['nr'] for a in itp_atoms if a['name'] in cap_all_names}
    cap_C_nr = next(a['nr'] for a in itp_atoms if a['name'] == cap_C_name)
    itp_root_nr = next(n for n in adj_nr[cap_C_nr] if n not in cap_nrs)

    # ── 3. BFS through base-RTP sidechain from CA ────────────────────────────
    BACKBONE = {'N', 'HN', 'H', 'CA', 'HA', 'C', 'O', 'OC1', 'OC2'}
    ffpath   = get_ff_path(ff)
    rtp_path = os.path.join(ffpath, 'merged.rtp')
    if not os.path.isfile(rtp_path):
        rtp_path = os.path.join(ffpath, 'aminoacids.rtp')
    rtp_entry     = RTPParser(rtp_path)[base_resname]
    rtp_atoms_dict = {a[0]: {'type': a[1]} for a in rtp_entry['atoms']}
    rtp_adj: Dict[str, set] = _dd(set)
    for bond in rtp_entry['bonds']:
        rtp_adj[bond[0]].add(bond[1])
        rtp_adj[bond[1]].add(bond[0])
    rtp_root = next((n for n in rtp_adj.get('CA', set()) if n not in BACKBONE), None)

    rtp_bfs: List[Tuple[str, str]] = []  # [(name, type)]
    if rtp_root:
        rtp_visited = {'CA'}
        rtp_q = deque([(rtp_root, 'CA')])
        while rtp_q:
            name, _ = rtp_q.popleft()
            if name in rtp_visited:
                continue
            rtp_visited.add(name)
            if name in rtp_atoms_dict:
                rtp_bfs.append((name, rtp_atoms_dict[name]['type']))
            for nb in sorted(rtp_adj.get(name, set())):
                if nb not in rtp_visited and nb in rtp_atoms_dict:
                    rtp_q.append((nb, name))

    # ── 4. Load sigma/eps data for type remapping ────────────────────────────
    def _is_abstract(t: str) -> bool:
        return t.startswith('AT_') or t.startswith('DUM_AT_')

    of_params: Dict[str, Tuple[float, float]] = {}
    if atomtypes_itp:
        for e in _parse_atomtypes_section(atomtypes_itp):
            of_params[e[0]] = (float(e[5]), float(e[6]))

    gaff_params: Dict[str, Tuple[float, float]] = {}
    for _cand_ff in ('gaff2', 'gaff'):
        try:
            _gpath = get_ff_path(_cand_ff)
            _gnb   = os.path.join(_gpath, 'ffnonbonded.itp')
            if os.path.isfile(_gnb):
                for e in _parse_atomtypes_section(_gnb):
                    gaff_params[e[0]] = (float(e[5]), float(e[6]))
                break
        except Exception:
            continue

    def _closest_gaff(sigma: float, eps: float) -> Optional[str]:
        if not gaff_params:
            return None
        return min(gaff_params,
                   key=lambda t: (gaff_params[t][0] - sigma)**2
                                + (gaff_params[t][1] - eps)**2)

    # ── 5. BFS through ITP sidechain, match against RTP BFS order ───────────
    def _elem(name: str) -> str:
        clean = name.lstrip('0123456789')
        alpha = ''.join(c for c in clean if c.isalpha())
        if not alpha:
            return name[0].upper()
        sym = alpha[0].upper()
        if len(alpha) > 1 and alpha[1].islower():
            sym += alpha[1]
        return sym

    itp_map: Dict[int, dict] = {}
    itp_visited = set(cap_nrs)
    itp_q: deque = deque([(itp_root_nr, cap_C_nr)])
    rtp_iter = iter(rtp_bfs)
    ligand_elem_count: Dict[str, int] = _dd(int)

    while itp_q:
        nr, _ = itp_q.popleft()
        if nr in itp_visited:
            continue
        itp_visited.add(nr)
        atom     = nr_to_atom[nr]
        raw_type = atom['type']
        rtp_match = next(rtp_iter, None)

        if rtp_match is not None:
            new_name, amber_type = rtp_match[0], rtp_match[1]
            gaff_type = (_AMBER_TO_GAFF.get(amber_type, raw_type)
                         if _is_abstract(raw_type) else raw_type)
        else:
            elem = _elem(atom['name'])
            ligand_elem_count[elem] += 1
            new_name   = f'{elem}{ligand_elem_count[elem]}'
            amber_type = None
            if _is_abstract(raw_type) and raw_type in of_params:
                s, e = of_params[raw_type]
                gaff_type = _closest_gaff(s, e) or raw_type
            else:
                gaff_type = raw_type

        itp_map[nr] = {'new_name': new_name, 'amber_type': amber_type,
                       'gaff_type': gaff_type}
        for nb in sorted(adj_nr.get(nr, set()),
                         key=lambda n: _elem(nr_to_atom[n]['name'])):
            if nb not in itp_visited:
                itp_q.append((nb, nr))

    # ── 6. Rewrite ITP with new names and types ──────────────────────────────
    old_to_new  = {nr_to_atom[nr]['name']: m['new_name']   for nr, m in itp_map.items()}
    old_to_gaff = {nr_to_atom[nr]['name']: m['gaff_type']  for nr, m in itp_map.items()}

    out_lines = []
    in_atoms  = False
    with open(itp_file) as fh:
        for raw in fh:
            stripped = raw.strip()
            if stripped.startswith('['):
                in_atoms = (stripped.strip('[]').strip() == 'atoms')
                out_lines.append(raw)
                continue
            if not in_atoms or not stripped or stripped.startswith(';'):
                out_lines.append(raw)
                continue
            parts = raw.rstrip('\n').split()
            if len(parts) < 7:
                out_lines.append(raw)
                continue
            try:
                int(parts[0])
            except ValueError:
                out_lines.append(raw)
                continue
            old_name = parts[4]
            parts[1] = old_to_gaff.get(old_name, parts[1])
            parts[4] = old_to_new.get(old_name, old_name)
            out_lines.append('   '.join(parts) + '\n')

    with open(out_itp, 'w') as fh:
        fh.writelines(out_lines)

    n_renamed = sum(1 for o, n in old_to_new.items()  if o != n)
    n_retyped = sum(1 for o, n in old_to_gaff.items() if o != n)
    logger.info("sanitize_itp: %d atoms renamed, %d types remapped → %s",
                n_renamed, n_retyped, out_itp)
    return out_itp, itp_map


def prepare_nsaa_ff(
    itp_file: str,
    base_resname: str,
    new_resname: str,
    ff: str,
    target_charge: int = 0,
    cap_C_override: Optional[str] = None,
    atomtypes_itp: Optional[str] = None,
    out_rtp: Optional[str] = None,
    out_supplement: Optional[str] = None,
) -> Tuple[str, str]:
    """Convert a CH3-capped ITP fragment into a proper force-field RTP entry.

    This is a **RDKit-free** alternative to :func:`mutate_nsaa` for the
    ITP-fragment workflow.  It derives the NSAA's residue topology directly
    from the ITP bond graph — no atom mapping heuristic is needed.

    Typical workflow
    ----------------
    1. Call this function once to register the NSAA with the force field::

           rtp, supp = prepare_nsaa_ff("sep.itp", "SER", "SEP", "amber14sbmut")

    2. Copy (or ``#include``) *supp* into the force-field directory so that
       ``pdb2gmx`` can resolve the cross-boundary AMBER↔GAFF bonded terms.

    3. Build the hybrid structure and topology with standard pmx tools::

           from pmx.alchemy import mutate
           mutate(m, mut_resid=42, mut_resname="SEP", ff="amber14sbmut",
                  inplace=True)
           m.write("mutant.pdb")
           # ... run pdb2gmx ...
           from pmx.forcefield import Topology
           from pmx.alchemy import gen_hybrid_top
           top = Topology("topol.top", ff="amber14sbmut", version="new")
           pmxtop, itps = gen_hybrid_top(top,
                                         extra_mtp_files=["S2SEP.mtp"],
                                         supplement_bonded_files=[supp])

    Parameters
    ----------
    itp_file : str
        GAFF/GAFF2/OpenFF ITP fragment with a CH3 cap at the CA position.
        For OpenFF or CGenFF inputs, call :func:`sanitize_itp` first (or
        pass *atomtypes_itp* to trigger auto-sanitization).
    base_resname : str
        Base amino acid name (e.g. ``'SER'``, ``'CYS'``).
    new_resname : str
        Name for the modified residue (e.g. ``'SEP'``, ``'CYSI'``).
    ff : str
        Force-field name.
    target_charge : int
        Desired integer total charge.  Any imbalance between the sum of
        backbone charges (from the base-AA RTP) plus sidechain charges
        (from the ITP) and *target_charge* is applied as a correction on CA.
    cap_C_override : str or None
        Override the auto-detected CH3 cap carbon name.
    atomtypes_itp : str or None
        Separate ``[ atomtypes ]`` ITP for OpenFF/CGenFF inputs.  When
        provided and the ITP contains abstract types (``AT_*``), auto-
        sanitization is run before processing.
    out_rtp : str or None
        Output RTP path.  Default: ``<ff_dir>/<new_resname>.rtp``.
    out_supplement : str or None
        Output bonded-supplement ITP path.
        Default: ``./<new_resname>_bonded.itp``.

    Returns
    -------
    out_rtp, out_supplement : str, str
        Paths to the written files.
    """
    from .ffparser import RTPParser

    ffpath = get_ff_path(ff)
    bCharmm = 'charmm' in ff.lower()

    if out_rtp is None:
        out_rtp = os.path.join(ffpath, f'{new_resname}.rtp')
    if out_supplement is None:
        out_supplement = f'./{new_resname}_bonded.itp'

    # ── 1. Base amino acid RTP ───────────────────────────────────────────────
    rtp_path = os.path.join(ffpath, 'merged.rtp')
    if not os.path.isfile(rtp_path):
        rtp_path = os.path.join(ffpath, 'aminoacids.rtp')
    rtp_parser  = RTPParser(rtp_path)
    base_entry  = rtp_parser[base_resname]
    base_atoms  = base_entry['atoms']   # [[name, type, charge, cgnr], ...]
    base_bonds  = base_entry['bonds']
    base_improps = base_entry.get('improps', [])

    bb_amber = {a[0]: {'type': a[1], 'charge': float(a[2])} for a in base_atoms}
    BACKBONE = {'N', 'HN', 'H', 'CA', 'HA', 'C', 'O', 'OT1', 'OT2', 'OC1', 'OC2'}

    # ── 2. Auto-sanitize ITP if abstract types detected ─────────────────────
    _raw = _parse_itp_atoms(itp_file)
    if (any(a['type'].startswith('AT_') for a in _raw) or
            len({a['name'] for a in _raw}) < len(_raw)):
        logger.info("prepare_nsaa_ff: sanitizing ITP (abstract/non-unique names)")
        itp_file, _ = sanitize_itp(itp_file, base_resname, ff,
                                   atomtypes_itp=atomtypes_itp,
                                   cap_C_override=cap_C_override)

    # ── 3. Parse ITP ────────────────────────────────────────────────────────
    itp_atoms     = _parse_itp_atoms(itp_file)
    _make_itp_names_unique(itp_atoms)
    itp_bonds_raw  = _parse_itp_bonds(itp_file)
    itp_bonds_full = _parse_itp_bonds_full(itp_file)
    itp_angles     = _parse_itp_angles(itp_file)
    itp_dihedrals  = _parse_itp_dihedrals(itp_file)

    nr_to_atom = {a['nr']: a for a in itp_atoms}
    nr_to_name = {a['nr']: a['name'] for a in itp_atoms}

    # ── 4. Identify CH3 cap ──────────────────────────────────────────────────
    cap_C_name, cap_H_names, connect_name, _ = _find_methyl_cap(
        itp_atoms, itp_bonds_raw)
    if cap_C_override is not None:
        cap_C_name = cap_C_override
        _adj_tmp = {a['nr']: set() for a in itp_atoms}
        for ai, aj in itp_bonds_raw:
            _adj_tmp[ai].add(aj); _adj_tmp[aj].add(ai)
        _cap_C_nr = next(a['nr'] for a in itp_atoms if a['name'] == cap_C_name)
        _h_nrs    = [n for n in _adj_tmp[_cap_C_nr]
                     if nr_to_atom[n]['name'].startswith('H')]
        cap_H_names  = [nr_to_atom[n]['name'] for n in _h_nrs]
        _non_h       = [n for n in _adj_tmp[_cap_C_nr] if n not in set(_h_nrs)]
        connect_name = nr_to_atom[_non_h[0]]['name'] if _non_h else connect_name

    if cap_C_name is None:
        raise ValueError(
            f"Cannot auto-detect CH3 cap in '{itp_file}'. "
            "Use cap_C_override to specify the cap carbon.")

    cap_all_names = {cap_C_name} | set(cap_H_names or [])
    cap_nrs       = {a['nr'] for a in itp_atoms if a['name'] in cap_all_names}

    # Sidechain atoms (ITP minus cap)
    sc_atoms  = [a for a in itp_atoms if a['nr'] not in cap_nrs]
    itp_type  = {a['name']: a['type'] for a in sc_atoms}   # name → GAFF type

    # ITP-internal adjacency (cap excluded) by atom name
    sc_adj: Dict[str, set] = {a['name']: set() for a in sc_atoms}
    for b in itp_bonds_full:
        ai, aj = b['indices']
        if ai in cap_nrs or aj in cap_nrs:
            continue
        na, nb = nr_to_name.get(ai), nr_to_name.get(aj)
        if na and nb:
            sc_adj[na].add(nb)
            sc_adj[nb].add(na)

    # Backbone adjacency (from base RTP bonds)
    bb_adj: Dict[str, set] = {a[0]: set() for a in base_atoms}
    for bond in base_bonds:
        bb_adj.setdefault(bond[0], set()).add(bond[1])
        bb_adj.setdefault(bond[1], set()).add(bond[0])

    bb_anchor = 'CA'   # cap_C stood in for CA; connect_name (e.g. CB) is its neighbour

    # AMBER type for connect_name in the base AA (needed for ffbonded lookups)
    base_connect_name = next(
        (n for n in bb_adj.get(bb_anchor, set())
         if n not in BACKBONE and n in bb_amber),
        None)
    connect_amber_type = (bb_amber[base_connect_name]['type']
                          if base_connect_name else itp_type.get(connect_name))
    connect_gaff_type  = itp_type.get(connect_name, connect_amber_type)

    # ── 5. Charge balance ────────────────────────────────────────────────────
    bb_q  = sum(float(a[2]) for a in base_atoms)
    sc_q  = sum(a['charge'] for a in sc_atoms)
    adj_q = target_charge - (bb_q + sc_q)
    logger.info(
        "prepare_nsaa_ff: backbone q=%.5f  sidechain q=%.5f  "
        "CA correction=%.5f  (target %d)",
        bb_q, sc_q, adj_q, target_charge)

    # ── 6. Write RTP ─────────────────────────────────────────────────────────
    header = _CHARMM_RTP_HEADER if bCharmm else _AMBER_RTP_HEADER
    max_bb_cgnr = max(int(a[3]) for a in base_atoms)
    rtp_lines = [header, f'\n[ {new_resname} ]\n\n', '[ atoms ]\n']

    for a in base_atoms:
        name, atype, charge, cgnr = a[0], a[1], float(a[2]), int(a[3])
        if name == bb_anchor:
            charge += adj_q
        rtp_lines.append(f'   {name:<6}   {atype:<15}  {charge:8.5f}  {cgnr}\n')
    for idx, a in enumerate(sc_atoms, start=1):
        rtp_lines.append(
            f'   {a["name"]:<6}   {a["type"]:<15}  {a["charge"]:8.5f}  '
            f'{max_bb_cgnr + idx}\n')

    rtp_lines.append('\n[ bonds ]\n')
    for bond in base_bonds:
        rtp_lines.append(f'   {bond[0]:<6}  {bond[1]:<6}\n')
    # cross-boundary bond (CA → connect_name / CB)
    if connect_name and connect_name not in bb_amber:
        rtp_lines.append(
            f'   {bb_anchor:<6}  {connect_name:<6}  '
            '; cross-boundary (cap–sidechain)\n')
    # ITP-internal bonds (no cap, no already-present backbone bonds)
    for b in itp_bonds_full:
        ai, aj = b['indices']
        if ai in cap_nrs or aj in cap_nrs:
            continue
        na, nb = nr_to_name.get(ai), nr_to_name.get(aj)
        if na and nb and not {na, nb} <= set(bb_amber):
            rtp_lines.append(f'   {na:<6}  {nb:<6}\n')

    if base_improps:
        rtp_lines.append('\n[ impropers ]\n')
        for imp in base_improps:
            rtp_lines.append('   ' + '  '.join(f'{x:<6}' for x in imp[:4]) + '\n')

    with open(out_rtp, 'w') as fh:
        fh.writelines(rtp_lines)
    logger.info("prepare_nsaa_ff: RTP written → %s", out_rtp)

    # ── 7. Write bonded supplement ITP ───────────────────────────────────────
    ffbonded = os.path.join(ffpath, 'ffbonded.itp')
    bond_db  = _parse_ffbonded_section(ffbonded, 'bondtypes',     2)
    angle_db = _parse_ffbonded_section(ffbonded, 'angletypes',    3)
    dihed_db = _parse_ffbonded_section(ffbonded, 'dihedraltypes', 4)

    def _lookup_bond(t1, t2):
        for e in bond_db:
            et = tuple(e['types'])
            if et in ((t1, t2), (t2, t1)):
                yield e

    def _lookup_angle(t1, t2, t3):
        for e in angle_db:
            et = tuple(e['types'])
            if et in ((t1, t2, t3), (t3, t2, t1)):
                yield e

    def _lookup_dihed(t1, t2, t3, t4):
        for e in dihed_db:
            et = tuple(e['types'])
            if et in ((t1, t2, t3, t4), (t4, t3, t2, t1)):
                yield e

    def _atype(name):
        """A-state (AMBER) type for any atom name."""
        return bb_amber.get(name, {}).get('type') or itp_type.get(name)

    def _btype(name):
        """B-state type: AMBER for backbone, GAFF for sidechain."""
        if name in bb_amber:
            return bb_amber[name]['type']
        return itp_type.get(name)

    seen_b, seen_a, seen_d = set(), set(), set()
    new_bonds, new_angles, new_diheds = [], [], []

    def _emit_bond(wt1, wt2, lt1, lt2):
        key = tuple(sorted([wt1, wt2]))
        if key in seen_b:
            return
        for e in _lookup_bond(lt1, lt2):
            if not e['params']:
                continue
            seen_b.add(key)
            new_bonds.append(([wt1, wt2], e['funct'], e['params']))
            # DUM_ variants for hybrid dummy-atom endpoints
            for dum_t, other in ((f'DUM_{wt1}', wt2), (f'DUM_{wt2}', wt1)):
                dk = tuple(sorted([dum_t, other]))
                if dk not in seen_b:
                    seen_b.add(dk)
                    new_bonds.append(([dum_t, other], e['funct'], e['params']))

    def _emit_angle(wt1, wt2, wt3, lt1, lt2, lt3):
        key = min((wt1, wt2, wt3), (wt3, wt2, wt1))
        if key in seen_a:
            return
        for e in _lookup_angle(lt1, lt2, lt3):
            if not e['params']:
                continue
            seen_a.add(key)
            new_angles.append(([wt1, wt2, wt3], e['funct'], e['params']))
            for i, dum_t in enumerate([f'DUM_{wt1}', f'DUM_{wt3}']):
                other = (wt3, wt2) if i == 0 else (wt1, wt2)
                dk = min((dum_t,) + other, other + (dum_t,))
                if dk not in seen_a:
                    seen_a.add(dk)
                    trip = ([dum_t, wt2, wt3] if i == 0
                            else [wt1, wt2, dum_t])
                    new_angles.append((trip, e['funct'], e['params']))

    def _emit_dihed(wt1, wt2, wt3, wt4, lt1, lt2, lt3, lt4):
        for e in _lookup_dihed(lt1, lt2, lt3, lt4):
            if not e['params']:
                continue
            funct = 9 if e['funct'] in (1, 9) else e['funct']
            mult  = (e['params'][2] if funct in (1, 9) and len(e['params']) >= 3
                     else '')
            key = (min((wt1, wt2, wt3, wt4), (wt4, wt3, wt2, wt1)),
                   funct, mult)
            if key not in seen_d:
                seen_d.add(key)
                new_diheds.append(([wt1, wt2, wt3, wt4], funct, e['params']))

    # ── 7a. Cross-boundary AMBER↔GAFF terms ──────────────────────────────────
    if connect_name and _atype(bb_anchor) and connect_amber_type:
        ca_at  = _atype(bb_anchor)     # AMBER type of CA
        cb_at  = connect_amber_type    # AMBER type of CB (for lookup)
        cb_bt  = connect_gaff_type     # GAFF  type of CB (for writing)

        _emit_bond(ca_at, cb_bt, ca_at, cb_at)

        for x in bb_adj.get(bb_anchor, set()):
            if _atype(x):
                _emit_angle(_atype(x), ca_at, cb_bt, _atype(x), ca_at, cb_at)

        for y in sc_adj.get(connect_name, set()):
            y_bt = _btype(y)
            y_at = itp_type.get(y, y_bt)  # use GAFF as best AMBER proxy
            if y_bt:
                _emit_angle(ca_at, cb_bt, y_bt, ca_at, cb_at, y_at)

        for x in bb_adj.get(bb_anchor, set()):
            x_at = _atype(x)
            for w in bb_adj.get(x, set()) - {bb_anchor}:
                w_at = _atype(w)
                if w_at and x_at:
                    _emit_dihed(w_at, x_at, ca_at, cb_bt,
                                w_at, x_at, ca_at, cb_at)

        for x in bb_adj.get(bb_anchor, set()):
            x_at = _atype(x)
            for y in sc_adj.get(connect_name, set()):
                y_bt = _btype(y)
                y_at = itp_type.get(y, y_bt)
                if x_at and y_bt:
                    _emit_dihed(x_at, ca_at, cb_bt, y_bt,
                                x_at, ca_at, cb_at, y_at)

    # ── 7b. ITP-internal GAFF terms (needed for pdb2gmx) ────────────────────
    nr_to_type = {a['nr']: a['type'] for a in itp_atoms}

    def _itp_types(indices):
        """Atom types for ITP term; None if any index is a cap atom."""
        types = []
        for nr in indices:
            if nr in cap_nrs:
                return None
            t = nr_to_type.get(nr)
            if t is None:
                return None
            types.append(t)
        return types

    for b in itp_bonds_full:
        types = _itp_types(b['indices'])
        if types is None or not b['params']:
            continue
        t1, t2 = types
        key = tuple(sorted([t1, t2]))
        if key in seen_b:
            continue
        seen_b.add(key)
        new_bonds.append(([t1, t2], b['funct'], b['params']))
        for dum_t, other in ((f'DUM_{t1}', t2), (f'DUM_{t2}', t1)):
            dk = tuple(sorted([dum_t, other]))
            if dk not in seen_b:
                seen_b.add(dk)
                new_bonds.append(([dum_t, other], b['funct'], b['params']))

    for ang in itp_angles:
        types = _itp_types(ang['indices'])
        if types is None or not ang['params']:
            continue
        t1, t2, t3 = types
        key = min((t1, t2, t3), (t3, t2, t1))
        if key in seen_a:
            continue
        seen_a.add(key)
        new_angles.append(([t1, t2, t3], ang['funct'], ang['params']))
        for dum_t, trip in ((f'DUM_{t1}', [f'DUM_{t1}', t2, t3]),
                            (f'DUM_{t3}', [t1, t2, f'DUM_{t3}'])):
            dk = min(tuple(trip), tuple(reversed(trip)))
            if dk not in seen_a:
                seen_a.add(dk)
                new_angles.append((trip, ang['funct'], ang['params']))

    for dih in itp_dihedrals:
        if dih['funct'] in (4, 2):   # improper — skip, handled by RTP
            continue
        types = _itp_types(dih['indices'])
        if types is None or not dih['params']:
            continue
        t1, t2, t3, t4 = types
        funct = 9 if dih['funct'] in (1, 9) else dih['funct']
        mult  = (dih['params'][2] if funct in (1, 9) and len(dih['params']) >= 3
                 else '')
        key = (min((t1, t2, t3, t4), (t4, t3, t2, t1)), funct, mult)
        if key not in seen_d:
            seen_d.add(key)
            new_diheds.append(([t1, t2, t3, t4], funct, dih['params']))

    # ── 7c. Optionally write DUM_* atomtypes for non-GAFF FFs ───────────────
    supp_lines = [
        f'; Supplementary bonded parameters: {base_resname} → {new_resname}\n',
        f'; Generated by pmx.nsaa.prepare_nsaa_ff — DO NOT EDIT MANUALLY\n\n',
    ]
    if atomtypes_itp is not None:
        at_entries = _parse_atomtypes_section(atomtypes_itp)
        if at_entries:
            supp_lines.append('[ atomtypes ]\n')
            supp_lines.append(
                '; name         at.num  mass         charge  ptype'
                '   sigma           epsilon\n')
            for parts in at_entries:
                name, at_num, mass, charge, ptype, sigma, epsilon = parts[:7]
                dum = f'DUM_{name}'
                supp_lines.append(
                    f'   {dum:<14}  {at_num:<4}  {mass:<12}  '
                    f'{charge:<7}  {ptype}  {sigma:<14}  {epsilon}\n')
            supp_lines.append('\n')

    if new_bonds:
        supp_lines.append('[ bondtypes ]\n')
        supp_lines.append('; type1    type2    funct  params\n')
        for types, funct, params in new_bonds:
            t1, t2 = types
            supp_lines.append(
                f'   {t1:<10}  {t2:<10}  {funct}   '
                + '   '.join(params) + '\n')
        supp_lines.append('\n')

    if new_angles:
        supp_lines.append('[ angletypes ]\n')
        supp_lines.append('; type1    type2    type3    funct  params\n')
        for types, funct, params in new_angles:
            t1, t2, t3 = types
            supp_lines.append(
                f'   {t1:<10}  {t2:<10}  {t3:<10}  {funct}   '
                + '   '.join(params) + '\n')
        supp_lines.append('\n')

    if new_diheds:
        supp_lines.append('[ dihedraltypes ]\n')
        supp_lines.append('; type1    type2    type3    type4    funct  params\n')
        for types, funct, params in new_diheds:
            t1, t2, t3, t4 = types
            supp_lines.append(
                f'   {t1:<10}  {t2:<10}  {t3:<10}  {t4:<10}  {funct}   '
                + '   '.join(params) + '\n')
        supp_lines.append('\n')

    with open(out_supplement, 'w') as fh:
        fh.writelines(supp_lines)

    logger.info("prepare_nsaa_ff: supplement written → %s", out_supplement)
    logger.info("  bonds: %d  angles: %d  dihedrals: %d",
                len(new_bonds), len(new_angles), len(new_diheds))
    logger.info(
        "  Next: copy '%s' into the FF dir so pdb2gmx can find the bonded params.",
        os.path.basename(out_supplement))

    return out_rtp, out_supplement
