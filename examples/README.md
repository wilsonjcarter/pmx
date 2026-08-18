# pmx Tutorials



Three self-contained examples covering the main pmx mutation workflows.
Each directory contains a `run.sh` (complete pipeline) and a `README.md`
(unique highlights + detailed walkthrough).

<table bgcolor="#ffffff" cellspacing="0" cellpadding="8" border="0">
<tr bgcolor="#ffffff">
<td align="center" width="33%" bgcolor="#ffffff">
<img src="imgs/schematic_p1.png" width="100%"/><br/>
<sub><b>01</b> &nbsp; amino acid substitution (W→F)</sub>
</td>
<td align="center" width="33%" bgcolor="#ffffff">
<img src="imgs/schematic_p2.png" width="100%"/><br/>
<sub><b>02</b> &nbsp; protonation state / pKa shift</sub>
</td>
<td align="center" width="33%" bgcolor="#ffffff">
<img src="imgs/schematic_p3.png" width="100%"/><br/>
<sub><b>03</b> &nbsp; disulfide bond formation</sub>
</td>
</tr>
</table>

## Examples

| # | Directory | System | What is calculated | FF | Charge change |
|---|-----------|--------|-------------------|----|---------------|
| 1 | `01_trpcage_W6F` | Trp Cage (1L2Y) | ΔΔG of W→F substitution | charmm36m-mut | none |
| 2 | `02_protonation_C2CM` | Thioredoxin (1ERT) | ΔpKa of Cys32 | charmm36m-mut | −1 (**doublebox**) |
| 3 | `03_disulfide_C2CD` | Thioredoxin (1ERT) | ΔE of disulfide formation | charmm36m-mut | none |

## Prerequisites

```bash
conda activate devpmx
source /usr/local/gromacs/bin/GMXRC   # or equivalent GROMACS activation
```

## Quick start

```bash
cd examples/01_trpcage_W6F
bash run.sh
```

## General pmx protein FEP workflow

```
1. pdb2gmx (wildtype)         → normalised atom names (.gro)
2. pmx mutate                 → hybrid structure (.pdb)
3. pdb2gmx (mutant)           → standard topology (.top)
4. pmx gentop                 → hybrid topology (pmxtop.top)
5a. gmx editconf              → simulation box (.gro)  [single-system only]
5b. gmx solvate + gmx genion  → solvated, ion-neutralised system
6.  endpoint MD + NEQ         → transition trajectories
7.  pmx analyze               → ΔG via Crooks/BAR
```

Step 5a (`gmx editconf`) defines the box for single-system setups (examples 1, 3). The
charge-changing example (2) uses `pmx doublebox` instead, which sets the box vectors
internally — `editconf` is not needed when using doublebox.


### Charge-changing mutations (doublebox)

When the alchemical mutation shifts the net charge of the system (example 2), a simple single-box simulation introduces finite-size PBC artefacts. The correct approach is the **single-box double-system** method:

```
protein leg:   CYS ──────────► CYM    Δq = −1  ┐
reference leg: CYM ──────────► CYS    Δq = +1  ┘ net box Δq = 0  ✓
```

The reference is typically the **biological peptide** in the same sequence context —
not a minimal capped tripeptide (except when this is the experimental reference). Using the real peptide captures
the free energy in the same sequence context and gives ΔΔG of binding directly. Use
`pmx doublebox` to place both structures in one box before solvation:

```bash
pmx doublebox -f1 processed.gro -f2 ref_processed.gro \
              -o doublebox.gro -r 2.5 -d 1.5
```

See the README in example 02 (§Step 5 or §Step 6) for the complete
reference-leg preparation and charge-neutral solvation.
