.. _ex_protonation:

Cysteine protonation state — thioredoxin Cys32 (1ERT)
======================================================

The active-site Cys32 of *E. coli* thioredoxin has a depressed thiol
pK\ :sub:`a` of ~6.3, compared to ~8.3 for a free cysteine in solution.
The electrostatic environment of the protein stabilises the thiolate (CYM)
relative to the protonated form (CYS).  This shift can be recovered from a
thermodynamic cycle involving alchemical deprotonation in both environments.

Thermodynamic cycle
-------------------

.. image:: cycles/cycle_protonation.png
    :width: 500px
    :align: center

.. math::

    \Delta\Delta G_\mathrm{prot} = \Delta G_\mathrm{CYS \to CYM}^\mathrm{protein}
                                  - \Delta G_\mathrm{CYS \to CYM}^\mathrm{water}

The pK\ :sub:`a` shift follows directly:

.. math::

    \Delta\mathrm{p}K_a = \frac{\Delta\Delta G_\mathrm{prot}}{RT \ln 10}

A negative :math:`\Delta\Delta G_\mathrm{prot}` means the protein stabilises
the thiolate — the pK\ :sub:`a` is depressed relative to solution.

Hybrid residue: ``C2CM`` (Cys → Cym).

- State A: CYS, SG type ``S``, charge −0.230, HG1 present
- State B: CYM, SG type ``SS``, charge −0.800, HG1 → dummy

Input files
-----------

:download:`1ERT.pdb <../../examples/02_protonation_C2CM/input/1ERT.pdb>` —
thioredoxin, reduced form, chain A (108 residues).

Setup
-----

::

    pmx mutate -f 1ERT.pdb -o mutant.pdb \
        --resid 32 --resname CM -ff charmm36m-mut

    gmx pdb2gmx -f mutant.pdb -o processed.gro -p topol.top \
        -ff charmm36m-mut -water tip3p -ignh -ss no

    pmx gentop -p topol.top -o pmxtop.top

For the reference leg, run the same transformation on a capped Cys dipeptide
(ACE-C-NME) in water.

Notes
-----

- ``-ss no`` prevents ``pdb2gmx`` from automatically forming a disulfide bond
  between Cys32 and Cys35; both remain free thiols in the reduced input structure.
- State B carries a net charge of −1 relative to state A.  The solution-leg
  simulation is straightforward to set up with a single counter-ion for the
  B state.  For the protein leg, carry a K\ :sup:`+` counter-ion that is
  present only in the B-state half of the free energy calculation.
- Cys35 (the other active-site cysteine) has a near-normal pK\ :sub:`a` and
  serves as an interesting comparison.
