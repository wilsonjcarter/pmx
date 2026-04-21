# pmx Example Cases

Self-contained example directories demonstrating different pmx mutation workflows.
Each subdirectory is independent and contains a `run.sh` script with the complete
workflow and a `README.md` explaining the biology and technical choices.

## Examples

| # | Directory | System | Mutation type | Force field |
|---|-----------|--------|---------------|-------------|
| 1 | `01_trpcage_W6A` | Trp Cage (1L2Y) | Standard AA substitution W→A | amber99sb-star-ildn-mut |
| 2 | `02_protonation_C2CM` | Thioredoxin (1ERT) | Cysteine protonation CYS→CYM | amber99sb-star-ildn-mut |
| 3 | `03_disulfide_C2CD` | Thioredoxin (1ERT) | Disulfide bond formation | amber99sb-star-ildn-mut |
| 4 | `04_nterm_deletion` | HSP90 MEEVD peptide (1ELR) | N-terminal residue deletion | charmm36m-mut |
| 5 | `05_phosphorylation_SEP` | SH2 domain substrate | NSAA: Ser→phosphoSer (SEP) | charmm36m-mut |

## Prerequisites

```bash
# GROMACS (gmx) must be in PATH
# pmx installed (pip install -e .) or conda
conda activate devpmx
source /usr/local/gromacs/bin/GMXRC

# For example 5 (NSAA), rdkit is NOT required if using prepare_nsaa_ff directly;
# but if you supply a custom ITP, have it parameterised beforehand.
```

## Quick start

```bash
cd examples/01_trpcage_W6A
bash run.sh
```

## General pmx protein FEP workflow

```
1. pmx mutate      → hybrid structure (.pdb)
2. gmx pdb2gmx     → standard topology (.top, .itp)
3. pmx gentop      → hybrid topology (pmxtop.top)
4. solvate + em/eq → production-ready system
5. pmx analyse     → ΔΔG from dH/dl or work files
```

For NSAA (example 5), step 1 is replaced by `pmx mutate_nsaa` which also
writes force-field patches (`.rtp`, `.mtp`, `_bonded.itp`) used in steps 2–3.
