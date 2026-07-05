# Disulfide H-torsion test: do the omitted torsions change ΔG?

The pre-graph disulfide hand-list (`_add_ss_bonded_terms` / `disulfide_morph`)
added the `SG-SG` bond, the `CB-SG-SG` angles and the `CA-CB-SG-SG'` dihedrals —
but **omitted the `HB-CB-SG-SG'` torsions** that a real GROMACS disulfide (and
the Level-2 graph rule, `induced_bond_morph`) include. Those torsions carry a
real CHARMM force constant (`SM SM CT2 HA2`, kφ = 0.661 kJ/mol). This example
measures their free-energy contribution directly.

Model: dimethyl disulfide (`CH3-S-S-CH3`), vacuum, ABCG2/GAFF2 params. The
`S-S` bond + its induced angles/dihedrals morph `k=0` (state A) → full (state B),
built two ways:

- **full** — the complete graph set from `induced_bond_morph` (7 dihedrals:
  `C-S-S-C` + 6 `H-C-S-S`);
- **heavy_only** — heavy-atom terms only (1 dihedral: `C-S-S-C`), reproducing
  the old hand-list.

We run reverse non-equilibrium switches for each and compare ΔG. The difference
is the free energy the hand-list omitted.

## Result (illustrative, short test run)

```
  full      : dG(A->B) = 21.0 kJ/mol
  heavy_only: dG(A->B) = 16.0 kJ/mol
  ==> H-torsion contribution  ddG = +5.0 kJ/mol = +1.2 kcal/mol
```

~1 kcal/mol — above the ~0.1–0.2 kcal/mol precision target of a well-converged
FEP, so the omission was a genuine systematic error, not noise.

## Run

```bash
export PMX_SRC=/path/to/pmx/src
bash run.sh                                   # defaults: NREP=10, 50 ps switches
NREP=100 NSTEPS_SW=200000 bash run.sh         # converged
```

Requirements: a single env with `gmx`, `antechamber`/`parmchk2`/`tleap`/`parmed`
(AmberTools ≥ 23 for ABCG2), `obabel`, and `pmx` importable.

## Caveats (important)

- **This is a vacuum DMDS model, not the shipped C2CD.** DMDS drops 6 H-torsions;
  a real protein `C2CD` disulfide drops only **4** (`HB-CB-SG-SG'`), and the
  protein backbone/geometry constrains the SG-SG torsions differently — so the
  real-system ΔΔG will differ (likely smaller). For the exact `C2CD` number,
  run the graph-vs-detector comparison through the full pmx CLI
  (`pdb2gmx → gen_hybrid_top`) on a 2-Cys protein — the deeper validation this
  toy motivates.
- Single-direction Jarzynski from few/short switches carries statistical noise;
  scale `NREP`/`NSTEPS_SW` for a tight number.
- Both ΔG values are the *bonded* free energy of turning the S-S terms on (no
  charge/LJ or thiol-H morph), which is what isolates the H-torsion effect.

## Files

| file | role |
|---|---|
| `run.sh` | runs both variants, reports ΔΔG |
| `build_hybrid.py` | dual-state topology; `heavy_only` drops H-torsions |
| `analyze.py` | Jarzynski ΔG |
| `mdp/` | em / equil / switch templates |
