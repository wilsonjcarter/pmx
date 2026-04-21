.. _ex_nterm_deletion:

N-terminal deletion — HSP90 MEEVD peptide (1ELR)
=================================================

The C-terminal MEEVD motif of HSP90 is the recognition sequence for TPR
co-chaperones.  The N-terminal methionine (Met1) makes contacts with the Hop
TPR2A domain, but its individual contribution to binding affinity is not easily
accessible by experiment.  Alchemical deletion gives a direct estimate of
:math:`\Delta\Delta G_\mathrm{bind}`.

Thermodynamic cycle
-------------------

.. image:: cycles/cycle_nterm_deletion.png
    :width: 500px
    :align: center

.. math::

    \Delta\Delta G_\mathrm{bind} = \Delta G_\mathrm{del}^\mathrm{complex}
                                  - \Delta G_\mathrm{del}^\mathrm{free\ peptide}

State A is the full MEEVD pentapeptide; state B is the truncated EEVD
tetrapeptide with Met1 fully dummified (zero charge, decoupled from
non-bonded interactions).

Hybrid residue: ``EdeN`` placed on Glu2 — the residue that becomes the
new N-terminus in state B.

- State A: standard Glu2 backbone (internal residue)
- State B: Glu2 with N-terminal patch charges; Met1 atoms → DUM\_ types

Input files
-----------

:download:`1ELR_peptide.pdb <../../examples/04_nterm_deletion/input/1ELR_peptide.pdb>` —
MEEVD pentapeptide, chain B of 1ELR, renumbered 1–5.

:download:`1ELR_complex.pdb <../../examples/04_nterm_deletion/input/1ELR_complex.pdb>` —
full Hop TPR2A + MEEVD complex (for the bound leg).

Prerequisites
-------------

The ``XdeN`` hybrid residue entries are generated once per force field::

    python3 -m pmx.scripts.generate_term_deletion \
        --ffdir $(python3 -c "from pmx.utils import get_ff_path; print(get_ff_path('charmm36m-mut'))") \
        --outdir .

This writes ``mutres_term.rtp`` and ``mutres_term.mtp`` to the current
directory.  Copy the RTP into the force-field directory so ``pdb2gmx`` can
find it, and pass the MTP to ``pmx gentop`` via ``--extra_mtp``.

Setup
-----

The deletion is triggered by targeting Met1 directly with ``--resname DEL``.
pmx detects that Met1 is the N-terminus, adds dummy protons to the new
N-terminus (Glu2), and renames Glu2 to ``EdeN`` internally::

    pmx mutate -f 1ELR_peptide.pdb -o mutant.pdb \
        --resid 1 --resname DEL -ff charmm36m-mut

    gmx pdb2gmx -f mutant.pdb -o processed.gro -p topol.top \
        -ff charmm36m-mut -water tip3p -ignh

    pmx gentop -p topol.top -o pmxtop.top \
        --extra_mtp mutres_term.mtp

Repeat for the complex leg using ``1ELR_complex.pdb`` as the starting
structure (with the same mutation applied to chain B).

Notes
-----

- Target the residue to be deleted (Met1, ``--resid 1 --resname DEL``), not
  its neighbor.  pmx looks up the neighbor automatically and renames it.
- Met1 atoms remain physically present in the structure throughout; they
  simply lose all non-bonded interactions in state B.  This is the standard
  pmx dual-topology approach.
- The complex leg requires that both the protein and the peptide chain are
  included in the topology; pass the separate ITP files to ``pdb2gmx``
  and ensure ``pmx gentop`` processes them with ``--recursive``.
