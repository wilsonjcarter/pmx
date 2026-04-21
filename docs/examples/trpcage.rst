.. _ex_trpcage:

Protein stability — Trp Cage W6A (1L2Y)
========================================

Tryptophan 6 packs against the hydrophobic core of the Trp Cage miniprotein
and is a major contributor to its folding free energy.
Here we compute how the W6A mutation affects thermodynamic stability.

Thermodynamic cycle
-------------------

The mutation free energy :math:`\Delta\Delta G_\mathrm{fold}` is obtained by
closing the thermodynamic cycle over two alchemical legs: the mutation carried
out in the **folded** protein, and the same mutation in a short **unfolded**
peptide (or single amino acid in solution) as the reference state.

.. image:: cycles/cycle_trpcage.png
    :width: 500px
    :align: center

.. math::

    \Delta\Delta G_\mathrm{fold} = \Delta G_\mathrm{mut}^\mathrm{folded}
                                  - \Delta G_\mathrm{mut}^\mathrm{unfolded}

A positive :math:`\Delta\Delta G_\mathrm{fold}` indicates destabilisation.
For W6A this is expected to be roughly +2 kcal/mol.

Input files
-----------

:download:`1L2Y.pdb <../../examples/01_trpcage_W6A/input/1L2Y.pdb>` — Trp
Cage miniprotein, NMR model 1, chain A (20 residues).

Setup
-----

Build the hybrid structure::

    pmx mutate -f 1L2Y.pdb -o mutant.pdb \
        --resid 6 --resname A -ff charmm36m-mut

Generate the standard topology::

    gmx pdb2gmx -f mutant.pdb -o processed.gro -p topol.top \
        -ff charmm36m-mut -water tip3p -ignh

Fill B states::

    pmx gentop -p topol.top -o pmxtop.top

Repeat for the reference leg using a Trp-Ala-Trp tripeptide (or an isolated
Trp dipeptide capped with ACE/CT3).  The difference of the two
:math:`\Delta G` values gives :math:`\Delta\Delta G_\mathrm{fold}`.

Notes
-----

- The Trp Cage has a single chain with no disulfides or non-standard residues,
  making it the simplest possible benchmark for the pmx workflow.
- 1L2Y is an NMR ensemble; only model 1 is used.  The structure already has
  hydrogens — pass ``-ignh`` to ``pdb2gmx`` to let it re-add them consistently
  with the chosen force field.
- A reasonable estimate of the unfolded-state :math:`\Delta G` can be obtained
  from a Trp-Ala capped dipeptide (ACE-W-A-NME) rather than the full Trp Cage
  sequence.
