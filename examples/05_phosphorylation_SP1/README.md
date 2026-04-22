# Alchemical phosphorylation: Ser → phosphoSer (SP1) in thioredoxin

<p align="center">
  <img src="../imgs/schematic_p5.png" alt="phospho thermocycle" width="800"/>
</p>

## TLDR

This example computes the phosphorylation ΔΔG (Ser → monoanionic phosphoserine SP1) using CHARMM36m's built-in SP1 residue — no force-field patching or custom ITP is needed, making it the fastest route to a phosphorylation free energy. The hybrid `S2P1` morphs the serine hydroxyl into a phosphate group by promoting five dummy phosphate atoms to real ones and converting `HG1` to a non-interacting dummy.

| | |
|---|---|
| **Hybrid** | `S2P1` |
| **State A** | Ser: `OG` type `OH1`, `HG1` present; charge 0 |
| **State B** | SP1: `OG` → `ON2`, `HG1` → `DUM_H`, dummy `DP/DO1P/DO2P/DOT/DHT` → real; charge −1 |
| **Charge shift** | 0 → −1 (add one K⁺ for state B) |
| **Special steps** | None — SP1 is already in the FF; use `P2` / `SP2` for dianionic (−2) form |

```bash
# Unique mutation step — see run.sh for the full pipeline
printf "42 P1\n" | pmx mutate \
    -f  wt.gro \
    -o  mutant.pdb \
    -ff charmm36m-mut
```

---

## Detailed walkthrough

## Step 1 — Prepare the wildtype topology

Run `pdb2gmx` with `-ignh` to strip and consistently rebuild hydrogens. The output `.gro` normalises atom names for `pmx mutate`; `wt.top` is discarded.

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

---

## Step 2 — Build the hybrid structure

Mutation code `P1` selects the `S2P1` hybrid (monoanionic phosphoserine). `pmx mutate` inserts the `OG`/`HG1` type changes and the five dummy phosphate atoms; adjust the residue number to your target Ser in pmx-renumbered coordinates.

```bash
printf "42 P1\n" | pmx mutate \
    -f      wt.gro \
    -o      mutant.pdb \
    -ff     charmm36m-mut
```

---

## Step 3 — Build the hybrid topology

Do **not** pass `-ignh` — `pmx mutate` has already positioned all atoms. `pdb2gmx` reads the `S2P1` residue definition from `mutres.rtp` (part of `charmm36m-mut`) and generates the standard GROMACS topology.

```bash
gmx pdb2gmx \
    -f      mutant.pdb \
    -o      processed.gro \
    -p      topol.top \
    -ff     charmm36m-mut \
    -water  tip3p
```

---

## Step 4 — Fill B-state bonded terms

`pmx gentop` reads the `S2P1` entry from `mutres.mtp` and fills in all B-state atom types, charges, and bonded parameters. State B carries one extra negative charge, requiring one additional K⁺ counter-ion.

```bash
pmx gentop \
    -p  topol.top \
    -o  pmxtop.top \
    -ff charmm36m-mut
```

Expected output:
```
log_> Hybrid Residue -> 42 | S2P1
log_> Making bonds for state B -> ...
log_> Total charge of state A = 0
log_> Total charge of state B = -1
```

---

## Step 5 — Solvate and add ions

```bash
gmx solvate -cp processed.gro -cs spc216.gro -p pmxtop.top -o solvated.gro

gmx grompp -f mdp/em.mdp -c solvated.gro -r solvated.gro \
           -p pmxtop.top -o ions.tpr -maxwarn 1
printf '13\n' | gmx genion -s ions.tpr -pname K -nname CL \
           -neutral -o ions.gro -p pmxtop.top
```

---

## Step 6 — Energy minimise

```bash
gmx grompp -f mdp/em.mdp -c ions.gro -r ions.gro \
           -p pmxtop.top -o em.tpr -maxwarn 1
gmx mdrun -v -deffnm em -ntmpi 1
```

---

## Step 7 — Equilibrate and run endpoint simulations

Phosphorylation introduces large electrostatic changes; allow at least 10 ns equilibration per endpoint and monitor phosphate–protein contacts before harvesting transition frames.

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

---

## Step 8 — Non-equilibrium transition simulations

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

## Step 9 — Analyse

`results.txt` gives ΔG1 (Ser→SP1 free energy cost in the protein). Repeat Steps 1–9 with a short Ser-containing reference peptide in water to get ΔG2, then ΔΔG = ΔG1 − ΔG2.

```bash
pmx analyze \
    -fA transitionA/frame*/ti.xvg \
    -fB transitionB/frame*/ti.xvg
```

```
ΔΔG = ΔG1 − ΔG2
```
