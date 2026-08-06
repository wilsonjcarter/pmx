"""Normalized target-residue definitions for covalent mutations.

A :class:`TargetResidue` wraps a pmx :class:`~pmx.model.Model` holding a single
residue whose atoms carry names, element symbols, and 3-D coordinates -- the
minimal geometry that pmx's graph-MCS mapping pipeline
(:func:`pmx.covalent.map_residues`) needs.  Force-field parameters (atom type, charge, mass) are
*not* required here: for an RTP target they are assigned downstream from the
force-field RTP, exactly as for a user-supplied structure.

The headline front end is :meth:`TargetResidue.from_rtp`, which builds a 3-D
residue from an RTP entry alone (connectivity + atom types), so a covalent
modification can be requested by *name* with no external structure file.

RDKit is required for :meth:`from_rtp`.
"""

from __future__ import annotations

import os
from functools import lru_cache
from typing import Dict, List, Optional, Tuple

import numpy as np

from ..atom import Atom
from ..model import Model
from ..ffparser import RTPParser, load_ff_rtp
from ..utils import get_ff_path

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    _RDKIT_AVAILABLE = True
except ImportError:  # pragma: no cover - exercised only without RDKit
    _RDKIT_AVAILABLE = False


# ── Element inference ─────────────────────────────────────────────────────────
#
# An RTP entry gives force-field *atom types*, not chemical elements.  The most
# robust element source is the atom-type mass (from ``atomtypes.atp``): every
# force field lists a mass close to a standard isotope-averaged atomic weight,
# so nearest-mass matching recovers the element without per-force-field tables.
# A first-letter fallback covers types missing from ``atomtypes.atp``.

_STD_MASS_ELEMENT: List[Tuple[float, str]] = [
    (1.008, "H"), (12.011, "C"), (14.007, "N"), (15.999, "O"),
    (18.998, "F"), (22.990, "Na"), (24.305, "Mg"), (30.974, "P"),
    (32.06, "S"), (35.45, "Cl"), (39.098, "K"), (40.078, "Ca"),
    (55.845, "Fe"), (65.38, "Zn"), (79.904, "Br"), (126.904, "I"),
    (78.971, "Se"),
]

# Multi-letter element names whose leading letters would otherwise be misread
# by the single-letter fallback (e.g. "CL" -> C).
_TWO_LETTER = {"CL": "Cl", "BR": "Br", "NA": "Na", "MG": "Mg",
               "CA": "Ca", "FE": "Fe", "ZN": "Zn", "SE": "Se"}


def _load_atomtype_masses(ffpath: str) -> Dict[str, float]:
    """Parse ``atomtypes.atp`` -> {atom_type: mass}."""
    masses: Dict[str, float] = {}
    atp = os.path.join(ffpath, "atomtypes.atp")
    if not os.path.isfile(atp):
        return masses
    with open(atp) as fh:
        for line in fh:
            line = line.split(";", 1)[0].split()
            if len(line) >= 2:
                try:
                    masses[line[0]] = float(line[1])
                except ValueError:
                    continue
    return masses


def _mass_to_element(mass: float) -> str:
    """Nearest standard atomic mass -> element symbol."""
    best, best_d = "C", 1e18
    for m, el in _STD_MASS_ELEMENT:
        d = abs(m - mass)
        if d < best_d:
            best, best_d = el, d
    return best


def _element_for_type(atomtype: str, masses: Dict[str, float]) -> str:
    """Infer element for a force-field atom type, mass first then name."""
    mass = masses.get(atomtype)
    if mass is not None and mass > 0.0:
        return _mass_to_element(mass)
    key = atomtype.upper()
    for two, el in _TWO_LETTER.items():
        if key.startswith(two):
            return el
    # DUM_* dummy types encode the shadowed type after the underscore
    if key.startswith("DUM_"):
        return _element_for_type(atomtype[4:], masses)
    return atomtype[0].upper() if atomtype else "C"


# ── RTP -> RDKit geometry ─────────────────────────────────────────────────────

@lru_cache(maxsize=8)
def _load_rtp(ff: str) -> RTPParser:
    """Parse and cache a force field's residue database (read-only use).

    ``RTPParser`` is slow on large merged CHARMM RTP files (tens of thousands of
    lines), and a single session may build several targets from the same force
    field; caching the parse keeps repeated ``from_rtp`` calls cheap.  Callers
    must treat the returned parser as read-only.

    Delegates to :func:`pmx.ffparser.load_ff_rtp`, so residues split out of the
    vendored file into a topic RTP (``nos.rtp``, ``caps.rtp``, ...) resolve here
    exactly as they already do for pdb2gmx.
    """
    return load_ff_rtp(ff)


# Residue-name synonyms so a generic/target name resolves to the force field's
# actual RTP entry: histidine tautomers (amber HID/HIE/HIP vs charmm
# HSD/HSE/HSP) and common protonation-state variants.  Ordered by preference;
# the histidine default (HID/HSD, N-delta) matches the pmx mutres convention.
_RESNAME_SYNONYMS = {
    "HIS":  ["HID", "HSD", "HIE", "HSE", "HIP", "HSP"],
    "HISD": ["HID", "HSD"], "HISE": ["HIE", "HSE"], "HISH": ["HIP", "HSP"],
    "HID":  ["HID", "HSD"], "HIE": ["HIE", "HSE"], "HIP": ["HIP", "HSP"],
    "HSD":  ["HSD", "HID"], "HSE": ["HSE", "HIE"], "HSP": ["HSP", "HIP"],
    "LYN":  ["LYN", "LSN"], "ASH": ["ASH", "ASPP"], "GLH": ["GLH", "GLUP"],
    # Both directions: CHARMM ships the disulfide-bonded cysteine as CYS2, Amber
    # as CYX.  Recipes name one endpoint, so without the reverse mapping a
    # disulfide targeting CYS2 is unresolvable on every Amber force field.
    "CYX":  ["CYX", "CYS2"], "CYS2": ["CYS2", "CYX"], "CYM": ["CYM"],
}


def _resolve_resname(resname: str, rtp) -> Optional[str]:
    """Resolve *resname* to a key present in *rtp*, via synonyms if needed."""
    if resname in rtp:
        return resname
    for alt in _RESNAME_SYNONYMS.get(resname.upper(), []):
        if alt in rtp:
            return alt
    return None


def _rtp_atoms_bonds(resname: str, ff: str):
    """Return (ffpath, resolved_resname, atom_rows, intra_bonds) for an RTP residue.

    ``atom_rows`` is the RTP ``[atoms]`` list [name, type, charge, cgnr];
    ``intra_bonds`` is the list of within-residue (name1, name2) pairs, with
    inter-residue ``+``/``-`` bonds dropped.
    """
    ffpath = get_ff_path(ff)
    rtp = _load_rtp(ff)
    resolved = _resolve_resname(resname, rtp)
    if resolved is None:
        searched = ", ".join(os.path.basename(s) for s in rtp.sources)
        raise KeyError(
            f"Residue '{resname}' not found in {ffpath} (searched: {searched}). "
            "Provide a structure/ITP target, or add the RTP entry.")
    entry = rtp[resolved]
    atom_rows = entry["atoms"]
    intra_bonds = [
        (a1, a2) for a1, a2 in entry["bonds"]
        if not a1.startswith(("+", "-")) and not a2.startswith(("+", "-"))
    ]
    return ffpath, resolved, atom_rows, intra_bonds


def _terminal_patch(atom_rows, intra_bonds, term: str, ff: str = None,
                    resname: str = None):
    """Apply the force field's terminal patch to RTP ``atom_rows``/``bonds``.

    Reuses :func:`pmx.nsaa._rtp_terminal` (the same patch the NSAA emit path
    uses) so a terminal residue is built from one definition, then returns the
    patched ``(atom_rows, intra_bonds)`` for embedding.  Only the atom and bond
    lists are needed here (elements + connectivity); the terminal impropers are
    resolved later, from the RTP, by :func:`build_hybrid_residue`.
    """
    from ..nsaa import _rtp_terminal
    entr = {
        "atoms": [list(r) for r in atom_rows],
        "bonds": [list(b) for b in intra_bonds],
        "diheds": [], "improps": [], "cmap": [],
    }
    entr = _rtp_terminal(entr, term, ff=ff, resname=resname)
    seen, bonds = set(), []
    for a1, a2 in entr["bonds"]:
        if a1.startswith(("+", "-")) or a2.startswith(("+", "-")):
            continue
        key = frozenset((a1, a2))
        if key in seen:
            continue
        seen.add(key)
        bonds.append((a1, a2))
    return entr["atoms"], bonds


def _embed_coords(elements: List[str],
                  bonds: List[Tuple[int, int]],
                  seed: int = 0xC0DE) -> List[List[float]]:
    """Build a single-bond RDKit mol from elements+connectivity and embed 3-D.

    Bond-order perception is deliberately skipped: an RTP entry lists every
    hydrogen and the full connectivity, so heavy-atom valence is already pinned
    by graph degree.  ETKDG embeds sane geometry directly on the single-bond
    graph -- avoiding fragile aromatic/conjugation perception.  (Aromatic rings
    embed slightly puckered; this only matters for genuinely new atoms and
    relaxes out under energy minimisation.)
    """
    if not _RDKIT_AVAILABLE:
        raise RuntimeError(
            "TargetResidue.from_rtp requires RDKit. Install with:\n"
            "    conda install -c conda-forge rdkit")

    rw = Chem.RWMol()
    for el in elements:
        rw.AddAtom(Chem.Atom(el))
    for i, j in bonds:
        rw.AddBond(i, j, Chem.BondType.SINGLE)
    mol = rw.GetMol()

    # An RTP entry lists every hydrogen explicitly, so no implicit Hs should be
    # added.  Where the explicit degree exceeds the element's neutral valence
    # (e.g. Lys NZ = N with 4 single bonds), assign the formal charge that makes
    # the graph valid (NH3+ -> +1).  This handles charged groups without any
    # bond-order perception.
    pt = Chem.GetPeriodicTable()
    for atom in mol.GetAtoms():
        atom.SetNoImplicit(True)
        dv = pt.GetDefaultValence(atom.GetAtomicNum())
        deg = atom.GetDegree()  # all bonds are single, so degree == bond order sum
        if dv > 0 and deg > dv:
            atom.SetFormalCharge(deg - dv)

    # Sanitize everything except valence/property checks so single-bonded rings
    # (e.g. Tyr/Trp) pass; ring perception is still needed for ETKDG.
    Chem.SanitizeMol(
        mol,
        sanitizeOps=(Chem.SanitizeFlags.SANITIZE_ALL
                     ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES),
    )

    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) != 0:
            raise RuntimeError(
                "RDKit failed to embed a 3-D conformer for the target residue.")

    conf = mol.GetConformer()
    return [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]


# ── Public API ────────────────────────────────────────────────────────────────

class TargetResidue:
    """A covalent-mutation target as a single-residue pmx Model with geometry.

    Attributes
    ----------
    model : pmx.model.Model
        Model containing exactly one residue.
    residue : pmx.molecule.Molecule
        The target residue (``model.residues[0]``).
    resname : str
        Residue name.
    source : str
        Front end that produced it (``'rtp'``, ``'structure'``, ...).
    """

    def __init__(self, model: Model, resname: str, source: str):
        self.model = model
        self.residue = model.residues[0]
        self.resname = resname
        self.source = source

    def __repr__(self) -> str:
        return (f"<TargetResidue {self.resname} "
                f"({len(self.residue.atoms)} atoms, source={self.source})>")

    # -- front ends ------------------------------------------------------------

    @classmethod
    def from_rtp(cls, resname: str, ff: str, seed: int = 0xC0DE,
                 term: Optional[str] = None,
                 verbose: bool = False) -> "TargetResidue":
        """Build a 3-D target residue from an RTP entry alone.

        Reads the ``[atoms]`` (names, force-field types, charges) and
        ``[bonds]`` (connectivity) of *resname* in the *ff* RTP database,
        infers elements from atom-type masses, and generates 3-D coordinates by
        RDKit distance-geometry embedding of the connectivity graph.

        Parameters
        ----------
        resname : str
            Residue name present in the force field's RTP (e.g. ``'SEP'``).
        ff : str
            Force-field name (resolved via ``GMXLIB``/``get_ff_path``).
        seed : int
            RDKit embedding random seed (deterministic geometry).
        term : {'N', 'C'}, optional
            Build the residue in its terminal-patched form (adds the terminus
            atoms and retypes/recharges the backbone via the force field's
            terminal patch).  Used to express a terminal deletion's neighbour as
            an ordinary internal->terminal residue morph.
        verbose : bool
            Print a short summary of the built residue.

        Returns
        -------
        TargetResidue
        """
        ffpath, resname, atom_rows, intra_bonds = _rtp_atoms_bonds(resname, ff)
        if term:
            atom_rows, intra_bonds = _terminal_patch(atom_rows, intra_bonds, term,
                                                     ff=ff, resname=resname)
        masses = _load_atomtype_masses(ffpath)

        # RTPParser returns RTP-native names (the RTP files and the mutres
        # library share this convention), so names are used directly.
        names = [row[0] for row in atom_rows]
        types = [row[1] for row in atom_rows]
        name2idx = {n: i for i, n in enumerate(names)}
        elements = [_element_for_type(t, masses) for t in types]

        bond_idx = [
            (name2idx[a1], name2idx[a2])
            for a1, a2 in intra_bonds
            if a1 in name2idx and a2 in name2idx
        ]

        coords = _embed_coords(elements, bond_idx, seed=seed)

        atoms = []
        for i, (name, atype, row) in enumerate(zip(names, types, atom_rows)):
            a = Atom(name=name, resname=resname, x=list(coords[i]),
                     symbol=elements[i], id=i + 1)
            a.resnr = 1
            a.chain_id = "A"
            a.atomtype = atype
            try:
                a.q = float(row[2])
            except (IndexError, ValueError, TypeError):
                a.q = 0.0
            a.m = masses.get(atype, 0.0)
            atoms.append(a)

        # attach intra-residue connectivity so downstream graph-based mapping
        # (pmx.covalent.mapping) has the molecular topology available
        for i, j in bond_idx:
            atoms[i].bonds.append(atoms[j])
            atoms[j].bonds.append(atoms[i])

        model = cls._model_from_atoms(atoms, resname)

        if verbose:
            print(f"TargetResidue.from_rtp: {resname} / {ff} -> "
                  f"{len(atoms)} atoms, {len(bond_idx)} bonds")

        return cls(model, resname, "rtp")

    @classmethod
    def from_structure(cls, path: str,
                       resname: Optional[str] = None) -> "TargetResidue":
        """Build a target residue from a PDB/GRO structure file.

        Parity with the existing ``nonstandardPDB`` input path: the residue's
        geometry is taken directly from the file.
        """
        model = Model()
        model.read(path)
        res = model.residues[0]
        return cls(model, resname or res.resname, "structure")

    @classmethod
    def from_itp(cls, itp_file: str, resname: str, seed: int = 0xC0DE,
                 verbose: bool = False) -> "TargetResidue":
        """Build a 3-D target residue from a GROMACS ITP fragment alone.

        Reads the ITP's ``[atoms]`` (names, types, charges, masses) and
        ``[bonds]`` (connectivity), infers elements from the atom masses, and
        generates 3-D coordinates by RDKit distance-geometry embedding of the
        bond graph -- exactly as :meth:`from_rtp` does for RTP residues.  An ITP
        carries no coordinates, so this makes a separate structure file
        unnecessary: the geometry is built from the connectivity, and the
        backbone is grafted from the WT during hybrid construction (the fragment
        is CH3-capped at CA).

        Parameters
        ----------
        itp_file : str
            GROMACS ITP fragment parameterizing the target residue.
        resname : str
            Name to assign the built residue (the hybrid target name).
        """
        from ..nsaa import _parse_itp_atoms, _parse_itp_bonds
        atom_rows = _parse_itp_atoms(itp_file)
        bonds_raw = _parse_itp_bonds(itp_file)
        if not atom_rows:
            raise ValueError(f"no [ atoms ] parsed from ITP: {itp_file}")

        nr2idx = {a['nr']: i for i, a in enumerate(atom_rows)}
        names = [a['name'] for a in atom_rows]
        elements = [_mass_to_element(a['mass']) if a['mass'] > 0.0
                    else _element_for_type(a['type'], {}) for a in atom_rows]
        bond_idx = [(nr2idx[i], nr2idx[j]) for i, j in bonds_raw
                    if i in nr2idx and j in nr2idx]

        coords = _embed_coords(elements, bond_idx, seed=seed)

        atoms = []
        for i, row in enumerate(atom_rows):
            a = Atom(name=names[i], resname=resname, x=list(coords[i]),
                     symbol=elements[i], id=i + 1)
            a.resnr = 1
            a.chain_id = "A"
            a.atomtype = row['type']
            a.q = row['charge']
            a.m = row['mass']
            atoms.append(a)

        for i, j in bond_idx:
            atoms[i].bonds.append(atoms[j])
            atoms[j].bonds.append(atoms[i])

        model = cls._model_from_atoms(atoms, resname)
        if verbose:
            print(f"TargetResidue.from_itp: {resname} <- {itp_file} -> "
                  f"{len(atoms)} atoms, {len(bond_idx)} bonds")
        return cls(model, resname, "itp")

    # -- helpers ---------------------------------------------------------------

    @staticmethod
    def _model_from_atoms(atoms: List[Atom], resname: str) -> Model:
        """Assemble a single-residue Model from a flat atom list."""
        model = Model(atoms=atoms, renumber_residues=False)
        for res in model.residues:
            res.resname = resname
        return model

    def write(self, path: str) -> None:
        """Write the target residue to a PDB/GRO file."""
        self.model.write(path)


def complete_target_from_structure(supplied, resname, ff, term=None):
    """Complete a target residue read from a structure against the force-field RTP.

    If *supplied* is MISSING atoms the RTP defines (e.g. a heavy-atom-only crystal
    template with no hydrogens), rebuild the full atom set (+ hydrogens, + types)
    from the RTP via :meth:`TargetResidue.from_rtp` and overlay only the supplied
    HEAVY-atom coordinates (Kabsch fit: matched heavy atoms snapped to the supplied
    positions exactly, the rest follow the fit and are refined downstream by
    RDKit/pdb2gmx). Structures already containing every RTP atom, and residues absent
    from the RTP, are returned unchanged.

    This is the single home for the "structure = heavy-atom geometry, RTP = identity/
    types/hydrogens" rule, shared by the covalent and nsaa mutation orchestrators.
    """
    try:
        full = TargetResidue.from_rtp(resname, ff, term=term).residue
    except Exception:
        return supplied
    sup = {a.name: a for a in supplied.atoms}
    if {a.name for a in full.atoms} <= set(sup):
        return supplied                                       # already complete
    is_heavy = lambda a: not a.name[:1].startswith(("H", "D"))
    common = [a for a in full.atoms if a.name in sup and is_heavy(a)]
    if len(common) < 3:
        return supplied
    P = np.array([list(a.x) for a in common])
    Q = np.array([list(sup[a.name].x) for a in common])
    Pc, Qc = P.mean(0), Q.mean(0)
    U, S, Vt = np.linalg.svd((P - Pc).T @ (Q - Qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    for a in full.atoms:
        if a.name in sup:
            a.x = list(sup[a.name].x)
        else:
            a.x = list(R @ (np.array(list(a.x)) - Pc) + Qc)
    full.set_resname(resname)
    return full
