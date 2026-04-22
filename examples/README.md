# pmx Example Cases

Six self-contained examples covering the main pmx mutation workflows.
Each directory contains a `run.sh` (complete pipeline) and a `README.md`
(unique highlights + detailed walkthrough).

## Examples

| # | Directory | System | What is calculated | FF |
|---|-----------|--------|-------------------|----|
| 1 | `01_trpcage_W6F` | Trp Cage (1L2Y) | ΔΔG of W→F substitution | charmm36m-mut |
| 2 | `02_protonation_C2CM` | Thioredoxin (1ERT) | pKa shift of Cys32 | charmm36m-mut |
| 3 | `03_disulfide_C2CD` | Thioredoxin (1ERT) | ΔG of disulfide formation | charmm36m-mut |
| 4 | `04_cterm_deletion` | HSP90 MEEVD peptide (1ELR) | ΔG of C-terminal truncation | charmm36m-mut |
| 5 | `05_phosphorylation_pTyr` | Lck SH2 domain (1AOT) | ΔΔG of Tyr phosphorylation (TP1) | charmm36m-mut |
| 6 | `06_nsaa_denovo` **(WIP)** | Protein with target Ser | ΔΔG of Ser→SEP via custom ITP | charmm36m-mut |

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
1. pdb2gmx (wildtype)  → normalised atom names (.gro)
2. pmx mutate          → hybrid structure (.pdb)
3. pdb2gmx (mutant)    → standard topology (.top)
4. pmx gentop          → hybrid topology (pmxtop.top)
5. solvate + ions + em → production-ready system
6. endpoint MD + NEQ   → transition trajectories
7. pmx analyze         → ΔG via Crooks/BAR
```

Example 4 adds a one-time `generate_term_deletion` step before step 2.
Example 6 replaces step 2 with `pmx mutate_nsaa` and adds `prepare_nsaa_ff` before it.
