.. _ex_disulfide:

Alchemical disulfide formation — thioredoxin Cys32–Cys35 (1ERT)
================================================================

The Cys32–Cys35 disulfide in *E. coli* thioredoxin defines its redox activity.
Starting from the reduced (dithiol) crystal structure, we compute the free
energy of disulfide **formation** inside the folded protein relative to a
model peptide in solution.

Thermodynamic cycle
-------------------

.. image:: cycles/cycle_disulfide.png
    :width: 500px
    :align: center

.. math::

    \Delta\Delta G_\mathrm{SS} = \Delta G_\mathrm{oxid}^\mathrm{protein}
                                - \Delta G_\mathrm{oxid}^\mathrm{peptide}

Hybrid residue: ``C2CD`` applied to **both** partner cysteines.

- State A: CYS, SG type ``S``, charge −0.230; thiol hydrogen present
- State B: CYS2 (disulfide form), SG type ``SM``, charge −0.080; HG1 → dummy

The cross-residue SG32–SG35 bond and its flanking angle and dihedral terms are
added automatically by ``pmx gentop`` and are not stored in the per-residue MTP.

Input files
-----------

:download:`1ERT.pdb <../../examples/03_disulfide_C2CD/input/1ERT.pdb>` —
thioredoxin, **reduced** form, chain A.  No disulfide bond is present in this
structure (SG32–SG35 distance ~3.9 Å).

Setup
-----

Mutate both cysteines in sequence::

    echo "32 CD" > mut.txt && pmx mutate -f 1ERT.pdb -o step1.pdb \
        --script mut.txt -ff charmm36m-mut

    echo "35 CD" > mut.txt && pmx mutate -f step1.pdb -o mutant.pdb \
        --script mut.txt -ff charmm36m-mut

Generate topology without forming a standard disulfide::

    gmx pdb2gmx -f mutant.pdb -o processed.gro -p topol.top \
        -ff charmm36m-mut -water tip3p -ignh -ss no

Fill B states — gentop detects the C2CD pair and injects the SG–SG bond::

    pmx gentop -p topol.top -o pmxtop.top

For the reference leg, use a Cys-Cys dipeptide (ACE-C-C-NME) in water.

Notes
-----

- The one-letter target code for the disulfide hybrid is ``CD``; this maps to
  the ``C2CD`` hybrid residue in the force field.
- The SG–SG distance in the reduced input structure is ~3.9 Å, versus ~2.05 Å
  in a disulfide.  The B state therefore involves significant bond compression;
  allow generous equilibration before collecting statistics.
- For the reverse transformation (disulfide breaking, i.e. reduction), use
  ``D2DC`` on both partner cysteines.
- Charge is conserved across states: both CYS and CYS2 carry zero net charge
  in the charmm36m force field.  No counter-ions are needed on this basis alone,
  though the overall protein charge should still be neutralised as normal.
