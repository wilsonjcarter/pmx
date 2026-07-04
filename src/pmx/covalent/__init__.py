"""Generalized covalent-modification framework for pmx.

This package unifies pmx's covalent-mutation machinery behind a single
abstraction: a *target residue* (named atoms + 3-D coordinates + force-field
parameters + connectivity) that can be produced from several front ends and
fed through one shared map -> hybridize -> emit pipeline.

Front ends for defining the target (see :mod:`pmx.covalent.target`):

* ``TargetResidue.from_rtp``       -- name an existing RTP residue; geometry is
                                      built from the RTP connectivity with RDKit.
* ``TargetResidue.from_structure`` -- supply a PDB/GRO structure for the residue.
* ``TargetResidue.from_itp``       -- supply a GROMACS ITP fragment (planned).
* ``TargetResidue.from_library``   -- use a pmx built-in residue (planned).

The downstream mapping (MCS + alignment), dummy-atom generation, and hybrid
RTP/MTP emission currently live in :mod:`pmx.nsaa`; this package is being grown
additively and will progressively host those pieces without changing behaviour.
"""

from .target import TargetResidue
from .mapping import map_residues, ResidueMapping, attach_bonds_from_rtp

__all__ = ["TargetResidue", "map_residues", "ResidueMapping",
           "attach_bonds_from_rtp"]
