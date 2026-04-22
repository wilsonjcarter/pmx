# [WIP] De-novo NSAA: Ser → phosphoSer (SEP) via ITP→RTP patching

> **Work in progress.** The `prepare_nsaa_ff` and `mutate_nsaa` tools are
> functional but the end-to-end workflow is still being validated.  Use
> `examples/06_phosphorylation_SP1/` for a production-ready phosphorylation
> calculation using the CHARMM36m built-in SP1 residue.

---

## When to use this workflow vs example 06

| Scenario | Recommended |
|----------|-------------|
| Phosphorylation using CHARMM36m SP1/SP2 | **Example 06** (fast, no ITP needed) |
| Phosphorylation with custom ITP parameters | **Example 07** (this example) |
| Any other PTM or non-standard AA | **Example 07** |

---

```
State A — unphosphorylated             State B — phosphorylated

  Ser42–OH   (protein)                   Ser42–OPO₃²⁻  (protein)
       |                                        |
      ΔG1 (protein)                           ΔG1
       |                                        |
  Ser42–OH   (water)    ─── ΔG2 (water) ─── Ser42–OPO₃²⁻  (water)

  ΔΔG = ΔG1 − ΔG2
```

We compute the free energy difference between the unphosphorylated (Ser) and
phosphorylated (phosphoSer / SEP) forms of a target serine.  Repeating the
calculation in water gives ΔG2; the difference ΔΔG is the effect of the
protein environment on the phosphorylation affinity.

This example uses the **NSAA (non-standard amino acid) workflow**, which
handles post-translational modifications and other residues not directly
encoded in the force-field libraries.  No RDKit is required; the ITP-to-RTP
framework (`prepare_nsaa_ff`) handles parameterisation directly.

**State A** — unphosphorylated: `OG` type `OH1` (charge −0.66), `HG1`
present (type `H`, charge 0.43)  
**State B** — phosphorylated: `OG` becomes the bridging oxygen of the
phosphate group, new phosphate atoms added, net charge shift of −2

---

## Prerequisites

In addition to the protein PDB you need two files for the SEP fragment:

1. **`input/sep.pdb`** — SEP structure with a methyl cap at the Cα position
   (the methyl replaces the backbone; Cα is the "cap carbon").  
   Source: RCSB CCD entry SEP, or build with Avogadro/Maestro/Open Babel.

2. **`input/sep.itp`** — GROMACS ITP fragment for SEP, parameterised with
   GAFF2 (e.g., via Antechamber/ACPYPE) or CGenFF (via ParamChem).  
   The cap methyl must match the cap in `sep.pdb`.

See `input/README_sep_params.md` for step-by-step instructions on preparing
these files.

---

## Step 1 — Prepare the wildtype topology

```bash
eval "$(python3 -c 'from pmx.gmx import set_gmxlib; import os; set_gmxlib(); print("export GMXLIB="+os.environ["GMXLIB"])')"

gmx pdb2gmx \
    -f      input/protein.pdb \
    -o      wt.gro \
    -p      wt.top \
    -ff     charmm36m-mut \
    -water  tip3p \
    -ignh
```

`-ignh` strips all existing hydrogens so GROMACS rebuilds them consistently.
The wildtype `.gro` is used only to normalise atom names before the NSAA
mutate step; the topology (`wt.top`) is discarded.

---

## Step 2 — Prepare the NSAA force-field patch

```bash
python3 -m pmx.nsaa.prepare_nsaa_ff \
    --itp       input/sep.itp \
    --base      SER \
    --name      SEP \
    --ff        charmm36m-mut \
    --charge    -2 \
    --out-rtp   SEP.rtp \
    --out-supp  SEP_bonded.itp
```

`prepare_nsaa_ff`:
- Reads the methyl-capped ITP fragment (`sep.itp` / `sep.pdb`)
- Builds a proper RTP entry for SEP using the Ser backbone plus SEP sidechain
  parameters
- Writes `SEP.rtp` into the FF directory so `pdb2gmx` can find it
- Writes `SEP_bonded.itp` — cross-boundary CHARMM↔GAFF bonded terms and
  ITP-internal GAFF bonded parameters needed for the phosphate sidechain
  atoms in state B

---

## Step 3 — Build the hybrid structure

```bash
pmx mutate_nsaa \
    -f        wt.gro \
    -o        mutant.pdb \
    --resid   42 \
    --resname SEP \
    -ff       charmm36m-mut \
    --itp     input/sep.itp \
    --nsapdb  input/sep.pdb
```

`pmx mutate_nsaa` replaces Ser42 with the hybrid residue `S2SEP`:
- State A: Ser sidechain (−OH group) — standard CHARMM36 Ser atoms
- State B: SEP sidechain (−OPO₃²⁻ group) — dummy phosphate atoms with zero
  charge/LJ in state A, real phosphate parameters in state B

Adjust `--resid 42` to the pmx-renumbered position of your target Ser.
`pmx mutate_nsaa` also writes `S2SEP.mtp` (the hybrid MTP for `pmx gentop`).

---

## Step 4 — Build the hybrid topology

```bash
gmx pdb2gmx \
    -f      mutant.pdb \
    -o      processed.gro \
    -p      topol.top \
    -ff     charmm36m-mut \
    -water  tip3p
```

Do **not** pass `-ignh` — `pmx mutate_nsaa` has already positioned all atoms.
`pdb2gmx` reads the `SEP` residue definition from `SEP.rtp` (placed in the FF
directory in Step 2) and generates the standard GROMACS topology.

---

## Step 5 — Fill B-state bonded terms

```bash
pmx gentop \
    -p                  topol.top \
    -o                  pmxtop.top \
    -ff                 charmm36m-mut \
    --extra_mtp         S2SEP.mtp \
    --supplement_bonded SEP_bonded.itp
```

`pmx gentop` reads `S2SEP.mtp` to locate the `S2SEP` hybrid residue and
`SEP_bonded.itp` to supply the GAFF/cross-boundary bonded parameters for the
phosphate group atoms in state B.

Expected output:
```
log_> Hybrid Residue -> 42 | S2SEP
log_> Making bonds for state B -> ...
log_> Total charge of state A = 0
log_> Total charge of state B = -2
```

State B carries two extra negative charges.  Two additional K⁺ counter-ions
are required for the state B simulation leg.

---

## Step 6 — Solvate and add ions

```bash
gmx solvate -cp processed.gro -cs spc216.gro -p pmxtop.top -o solvated.gro

gmx grompp -f mdp/em.mdp -c solvated.gro -r solvated.gro \
           -p pmxtop.top -o ions.tpr -maxwarn 1
printf '13\n' | gmx genion -s ions.tpr -pname K -nname CL \
           -neutral -o ions.gro -p pmxtop.top
```

Neutralise relative to state A.  For state B you will need two additional K⁺
ions (or apply a charge-correction in post-processing).

---

## Step 7 — Energy minimise

```bash
gmx grompp -f mdp/em.mdp -c ions.gro -r ions.gro \
           -p pmxtop.top -o em.tpr -maxwarn 1
gmx mdrun -v -deffnm em -ntmpi 1
```

---

## Step 8 — Equilibrate and run endpoint simulations

```bash
mkdir -p stateA
gmx grompp -f mdp/eqA.mdp -c em.gro -r em.gro \
           -p pmxtop.top -o stateA/md.tpr -maxwarn 1
gmx mdrun -v -deffnm stateA/md

mkdir -p stateB
gmx grompp -f mdp/eqB.mdp -c em.gro -r em.gro \
           -p pmxtop.top -o stateB/md.tpr -maxwarn 1
gmx mdrun -v -deffnm stateB/md
```

Phosphorylation introduces large electrostatic changes; allow at least 10 ns
equilibration per endpoint.  Monitor the phosphate–protein contacts to confirm
the local environment has stabilised before harvesting transition frames.

---

## Step 9 — Non-equilibrium transition simulations

```bash
mkdir -p transitionA
printf '0\n' | gmx trjconv -f stateA/md.trr -s stateA/md.tpr \
    -b 5000 -sep -o transitionA/frame.pdb

for i in $(seq 0 99); do
    mkdir -p transitionA/frame${i}
    mv transitionA/frame${i}.pdb transitionA/frame${i}/initial.pdb
    gmx grompp -f mdp/tiA.mdp \
        -c transitionA/frame${i}/initial.pdb \
        -r transitionA/frame${i}/initial.pdb \
        -p pmxtop.top -o transitionA/frame${i}/ti.tpr -maxwarn 1
    gmx mdrun -v -deffnm transitionA/frame${i}/ti
done

mkdir -p transitionB
printf '0\n' | gmx trjconv -f stateB/md.trr -s stateB/md.tpr \
    -b 5000 -sep -o transitionB/frame.pdb

for i in $(seq 0 99); do
    mkdir -p transitionB/frame${i}
    mv transitionB/frame${i}.pdb transitionB/frame${i}/initial.pdb
    gmx grompp -f mdp/tiB.mdp \
        -c transitionB/frame${i}/initial.pdb \
        -r transitionB/frame${i}/initial.pdb \
        -p pmxtop.top -o transitionB/frame${i}/ti.tpr -maxwarn 1
    gmx mdrun -v -deffnm transitionB/frame${i}/ti
done
```

---

## Step 10 — Analyse

```bash
pmx analyze \
    -fA transitionA/frame*/ti.xvg \
    -fB transitionB/frame*/ti.xvg
```

`results.txt` gives ΔG1 (Ser→SEP free energy cost in the protein).  Repeat
Steps 1–10 with a short Ser-containing reference peptide in water to get ΔG2,
then:

```
ΔΔG = ΔG1 − ΔG2
```

A negative ΔΔG means phosphorylation is more favourable in the protein
environment than in water (i.e. the protein stabilises the phosphorylated
form, as seen in many regulatory phosphorylation sites).

---

## Notes

- The NSAA workflow is generic: replace `SER`/`SEP` with any canonical AA
  and its PTM analogue.  The `--base` argument tells `prepare_nsaa_ff` which
  canonical backbone to use.
- The charge difference between state A (0) and state B (−2) requires a
  charge-correction term for quantitative free energies in periodic boundary
  conditions.  See Reif & Oostenbrink (2014) for the analytical correction.
- `pmx mutate_nsaa` does not require RDKit.  If you want to generate the
  ITP from SMILES automatically, install RDKit and use `pmx prepare_nsaa`.
- The phosphate parameters in `sep.itp` should be compatible with the protein
  force field.  GAFF2 phosphate parameters work well with CHARMM36m for
  mixed-FF calculations; CGenFF parameters are another option for a fully
  consistent treatment.
- Adjust `--resid 42` to match the target serine in your own structure (pmx
  numbers residues from 1 within each chain by default).
- For the simpler case where the target residue IS already in the force-field
  library (SP1/SP2 for phosphoserine), see `examples/06_phosphorylation_SP1/`.
