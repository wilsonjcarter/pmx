.. _ex_phosphorylation:

Serine phosphorylation — SEP via NSAA workflow
===============================================

Phosphorylation of serine introduces a net −2 charge and dramatically reshapes
the electrostatic surface of a protein.  SH2 and 14-3-3 domains are known to
bind phosphopeptides with nanomolar affinity, whereas the unphosphorylated
sequence binds weakly if at all.  The free energy cost of phosphorylation —
and its coupling to binding — can be decomposed using the thermodynamic cycle
below.

Thermodynamic cycle
-------------------

.. image:: cycles/cycle_phosphorylation.png
    :width: 500px
    :align: center

.. math::

    \Delta\Delta G_\mathrm{bind} = \Delta G_\mathrm{phospho}^\mathrm{complex}
                                  - \Delta G_\mathrm{phospho}^\mathrm{free}

State A is the serine-containing peptide (neutral sidechain, −OH); state B
is phosphoserine (SEP, monoanionic −OPO\ :sub:`3`\ H\ :sup:`−` at the
parameters encoded here, net −1 from the force-field entry, or −2 for the
dianionic form depending on pH and parameterisation).

Hybrid residue: ``S2SEP`` (Ser → SEP).

This example uses the **NSAA workflow** because SEP is not a standard amino
acid and requires a force-field patch before ``pdb2gmx`` can process it.

Input files
-----------

:download:`protein.pdb <../../examples/05_phosphorylation_SEP/input/protein.pdb>` —
Src SH2 domain, chain A of 2O88 (451 atoms).

You will also need a GAFF2-parameterised ITP for SEP (``sep.itp``) and a
capped reference structure (``sep.pdb``).  See
:download:`README_sep_params.md <../../examples/05_phosphorylation_SEP/input/README_sep_params.md>`
for how to generate these with Antechamber/ACPYPE.

Setup
-----

**Step 1** — build the force-field RTP entry and bonded supplement from the
ITP fragment::

    # Python API
    from pmx.nsaa import prepare_nsaa_ff

    out_rtp, out_supplement = prepare_nsaa_ff(
        itp_file      = 'sep.itp',
        base_resname  = 'SER',
        new_resname   = 'SEP',
        ff            = 'charmm36m-mut',
        target_charge = -1,          # monoanionic at physiological pH
        out_rtp       = 'SEP.rtp',
        out_supplement= 'SEP_bonded.itp',
    )

**Step 2** — build the hybrid structure::

    pmx mutate_nsaa -f protein.pdb -o mutant.pdb \
        --resid 42 --resname SEP -ff charmm36m-mut \
        --itp sep.itp --nsapdb sep.pdb

**Step 3** — standard GROMACS topology::

    gmx pdb2gmx -f mutant.pdb -o processed.gro -p topol.top \
        -ff charmm36m-mut -water tip3p -ignh

**Step 4** — fill B states, injecting the cross-boundary AMBER↔GAFF bonded
terms from the supplement::

    pmx gentop -p topol.top -o pmxtop.top \
        --extra_mtp S2SEP.mtp \
        --supplement_bonded SEP_bonded.itp

Notes
-----

- ``prepare_nsaa_ff`` reads the CH\ :sub:`3`-capped ITP and maps sidechain
  atoms onto the canonical Ser backbone from the RTP.  The backbone atoms
  (N, HN, CA, HA, C, O) are taken verbatim from the Ser entry; everything
  beyond CB comes from the ITP.
- The bonded supplement (``SEP_bonded.itp``) provides bonded parameters at
  the CA–CB boundary where CHARMM (backbone) and GAFF (sidechain) atom types
  meet.  These cross-type terms are not in the standard CHARMM ``ffbonded.itp``
  and must be injected before the B-state parameter search runs.
- State B gains −1 or −2 charge depending on whether the monoanionic or
  dianionic form of the phosphate is used.  Add counter-ions accordingly for
  the B-state leg.
- The same workflow applies to phosphothreonine (TPO) and phosphotyrosine
  (PTR), which have complete RTP entries in ``charmm36m-mut`` and therefore
  require no ``prepare_nsaa_ff`` call — only the MTP and bonded supplement
  steps remain.
