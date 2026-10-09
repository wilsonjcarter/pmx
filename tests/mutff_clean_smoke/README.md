# mutff_clean smoke harness

Validation harness for the clean force-field set in `src/pmx/data/mutff_clean/`.
**Relocated here (2026-09-25)** from that data dir — it is tooling, not FF content, so it
does not belong shipped inside the packaged FF directory.

## What it does
For each PTM/covalent target, builds a capped peptide `ACE-ALA-<WT>-ALA-NME`, runs
`pmx mutate --generative` to the target, then `pdb2gmx -> pmx gentop -> grompp`, and reports
any missing FF parameters. Matrix per FF = phospho {SER,THR,TYR}×{p1,p2} + methyl LYS×{me1,me2,me3}
+ (amber) TYM + GACE + (charmm) TYM.

## Run
```bash
conda activate devpmx           # pmx + rdkit
source /usr/local/gromacs/bin/GMXRC
bash tests/mutff_clean_smoke/smoke_test.sh <ff-name>   # e.g. amber14sb-q-mut, charmm36m-ecc-mut
```
The script sets `GMXLIB` to `../../src/pmx/data/mutff_clean` itself and forces
`PYTHONPATH` to `Documents/pmx/src` (so it never picks up the off-limits GitHub/pmx install).
Build outputs land in `_smoke/` here (git-ignored).

## Current status
11 FFs, **119/119 grompp-clean** (9 amber×11 + 2 charmm×10). Full build log + provenance:
`src/pmx/data/mutff_clean/HANDOFF.md`.

Files: `smoke_test.sh` (driver), `_mk_wt_peptide.py` (builds the capped WT peptide),
`_mk_gace_pep.py` (builds the GSH `GACE-CYS-GLY-NME` backbone for the GACE case).
