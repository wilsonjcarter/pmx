"""Normalized target-residue definitions for covalent mutations.

A :class:`TargetResidue` wraps a pmx :class:`~pmx.model.Model` holding a single
residue whose atoms carry names, element symbols, and 3-D coordinates -- the
minimal geometry that pmx's mapping pipeline (``pmx.nsaa.domapping`` /
``mutate_nsaa``) needs.  Force-field parameters (atom type, charge, mass) are
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

from ..atom import Atom
from ..model import Model
from ..ffparser import RTPParser
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
def _load_rtp(rtp_file: str) -> RTPParser:
    """Parse and cache an RTP database (read-only use).

    ``RTPParser`` is slow on large merged CHARMM RTP files (tens of thousands of
    lines), and a single session may build several targets from the same force
    field; caching the parse keeps repeated ``from_rtp`` calls cheap.  Callers
    must treat the returned parser as read-only.
    """
    return RTPParser(rtp_file)


def _rtp_atoms_bonds(resname: str, ff: str):
    """Return (ffpath, atom_rows, intra_bonds) for an RTP residue.

    ``atom_rows`` is the RTP ``[atoms]`` list [name, type, charge, cgnr];
    ``intra_bonds`` is the list of within-residue (name1, name2) pairs, with
    inter-residue ``+``/``-`` bonds dropped.
    """
    ffpath = get_ff_path(ff)
    rtp_file = os.path.join(ffpath, "merged.rtp")
    if not os.path.isfile(rtp_file):
        rtp_file = os.path.join(ffpath, "aminoacids.rtp")
    rtp = _load_rtp(rtp_file)
    if resname not in rtp:
        raise KeyError(
            f"Residue '{resname}' not found in {rtp_file}. "
            "Provide a structure/ITP target, or add the RTP entry.")
    entry = rtp[resname]
    atom_rows = entry["atoms"]
    intra_bonds = [
        (a1, a2) for a1, a2 in entry["bonds"]
        if not a1.startswith(("+", "-")) and not a2.startswith(("+", "-"))
    ]
    return ffpath, atom_rows, intra_bonds


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
        verbose : bool
            Print a short summary of the built residue.

        Returns
        -------
        TargetResidue
        """
        ffpath, atom_rows, intra_bonds = _rtp_atoms_bonds(resname, ff)
        masses = _load_atomtype_masses(ffpath)

        # RTPParser applies library._aliases (e.g. HB1 -> 1HB) for standard
        # residues, but the RTP files (and the mutres library) use canonical
        # names, and residues without an alias entry (e.g. SEP) are left as-is.
        # Undo the alias so the engine emits consistent RTP-native names.
        from ..library import _aliases
        _unalias = {v: k for k, v in _aliases.get(resname, {}).items()}

        def _canon(n):
            return _unalias.get(n, n)

        names = [_canon(row[0]) for row in atom_rows]
        types = [row[1] for row in atom_rows]
        name2idx = {n: i for i, n in enumerate(names)}
        elements = [_element_for_type(t, masses) for t in types]

        bond_idx = [
            (name2idx[_canon(a1)], name2idx[_canon(a2)])
            for a1, a2 in intra_bonds
            if _canon(a1) in name2idx and _canon(a2) in name2idx
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
    def from_itp(cls, *args, **kwargs) -> "TargetResidue":
        """Planned: build a target residue from a capped ITP fragment.

        The capped-fragment logic currently lives in
        ``pmx.nsaa._prepare_nsaa_fragment`` / ``_find_methyl_cap``; it will be
        lifted here as the ITP front end.
        """
        raise NotImplementedError(
            "from_itp is not implemented yet; use pmx.nsaa.mutate_nsaa with an "
            "itp= argument for the capped-fragment workflow.")

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
