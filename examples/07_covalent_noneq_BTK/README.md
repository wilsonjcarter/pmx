# BTK-style covalent inhibitor: alchemical thioether formation (vacuum)

A minimal, runnable demonstration of the pmx covalent framework: form the
covalent **acrylamide → Cys thioether** bond alchemically and measure the
non-equilibrium switching free energy — with the alchemical bonded
transformation defined entirely by `pmx.covalent.edit.induced_bond_morph`
(the Level-2 graph rule) and ligand parameters from **AmberTools/antechamber
ABCG2 + GAFF2**.

Model system: `CH3-S-CH2-CH2-C(=O)-NH-CH3` (SMILES `CSCCC(=O)NC`), the Michael
adduct of a methanethiol Cys-SG mimic and an N-methyl acrylamide warhead. The
forming bond is `S–C1` (Cys sulfur → former acrylamide β-carbon).

## What it does

1. `obabel` → 3-D structure; `antechamber -c abcg2 -at gaff2` → GAFF2 params +
   ABCG2 charges; `parmchk2`/`tleap`/`parmed` → GROMACS state-B topology.
2. `build_hybrid.py` builds a **dual-state (A/B)** topology: the forming `S–C1`
   bond and *its graph-induced angles/dihedrals* (selected by
   `induced_bond_morph`) morph `k = 0` (state A) ↔ full GAFF value (state B).
   Charges/LJ are unperturbed, so `dH/dλ` is purely the covalent bonded morph.
3. Minimize + equilibrate the well-defined **state B** (bonded thioether).
4. Run `NREP` **reverse** non-equilibrium switches (λ: 1→0, bond turns off),
   collect the work of each.
5. `analyze.py` → Jarzynski ΔG.

We switch *reverse* (from the bonded, intact state B) on purpose: the state-A
end (bond off) would let the two fragments drift apart, so state B is the
well-defined ensemble to sample. This is exactly the problem the split-atom /
dummy scheme solves for a rigorous bidirectional estimate (see the disulfide
demo).

## Requirements

A single environment on your machine with **all** of:
- `gmx` (GROMACS 2020+)
- `antechamber`, `parmchk2`, `tleap` (AmberTools ≥ 23 for ABCG2), `parmed`
- `obabel`
- `pmx` importable (or set `PMX_SRC=/path/to/pmx/src`, with its `_pmx`
  extension built for that Python)

## Run

```bash
export PMX_SRC=/path/to/pmx/src          # if pmx isn't pip-installed
bash run.sh                              # defaults: NREP=10, 50 ps switches
# scale up for convergence:
NREP=50 NSTEPS_SW=200000 NSTEPS_EQ=2000000 bash run.sh
```

### Choosing the forming-bond atoms for a different ligand

`build_hybrid.py` takes the two forming-bond atoms **by name**. antechamber/
GROMACS names atoms by element in SMILES order (`C, S, C1, C2, …`). Inspect the
`[ atoms ]` block of `state_b.top` (or run `build_hybrid.py` and read its
summary) and set:

```bash
BOND_A=S BOND_B=C1 SMILES="CSCCC(=O)NC" bash run.sh
```

## Caveats (honest)

- Single-direction Jarzynski is **biased** for few/fast switches. For a
  converged ΔG, run both directions and use BAR/Crooks — which needs a
  well-defined state-A ensemble, i.e. the **split-atom/dummy** hybrid.
- This isolates the *bonded* free energy of bond formation (no charge/LJ or
  thiol-H/hybridization morph). It demonstrates the machinery and the graph
  rule producing a GROMACS-valid, minimizable, switchable covalent topology; it
  is not a full covalent-inhibitor binding free energy.
- Vacuum: a single molecule in a large PBC box with plain cutoffs.

## Files

| file | role |
|---|---|
| `run.sh` | end-to-end driver |
| `build_hybrid.py` | dual-state topology from `induced_bond_morph` |
| `analyze.py` | Jarzynski ΔG from the switching works |
| `mdp/` | em / equil / switch templates |
