# [WIP] De-novo NSAA: Ser → phosphoSer (SEP) via ITP→RTP patching

> **Work in progress.** The `prepare_nsaa_ff` and `mutate_nsaa` tools are
> functional but the end-to-end workflow is still being validated.  Use
> `examples/05_phosphorylation_SP1/` for a production-ready phosphorylation
> calculation using the CHARMM36m built-in SP1 residue.

## TLDR

This example computes the phosphorylation ΔΔG (Ser → dianionc phosphoserine SEP) using a custom ITP fragment rather than a built-in FF residue — the `prepare_nsaa_ff` ITP→RTP patching step is the unique feature and generalises to **any** non-standard amino acid or PTM. Unlike example 05, this workflow uses `pmx mutate_nsaa` and requires `--extra_mtp`/`--supplement_bonded` flags in `pmx gentop` to supply cross-boundary bonded terms.

| | |
|---|---|
| **Hybrid** | `S2SEP` |
| **State A** | Ser: `OG` type `OH1`, `HG1` present; charge 0 |
| **State B** | SEP: bridging phosphate oxygen, full phosphate group; charge −2 |
| **Charge shift** | 0 → −2 (add two K⁺ for state B) |
| **Special steps** | (a) `prepare_nsaa_ff` to patch FF from `sep.itp`; (b) `pmx mutate_nsaa` instead of `pmx mutate`; (c) `--extra_mtp S2SEP.mtp --supplement_bonded SEP_bonded.itp` in `pmx gentop` |

```bash
# Unique steps — see run.sh for the full pipeline

# 1. Patch the force field from the custom ITP
python3 -m pmx.nsaa.prepare_nsaa_ff \
    --itp sep.itp --base SER --name SEP \
    --ff charmm36m-mut --charge -2 \
    --out-rtp SEP.rtp --out-supp SEP_bonded.itp

# 2. Build the hybrid structure with mutate_nsaa
pmx mutate_nsaa \
    -f wt.gro -o mutant.pdb \
    --resid 42 --resname SEP \
    -ff charmm36m-mut \
    --itp input/sep.itp --nsapdb input/sep.pdb
```

---

## When to use this workflow vs example 05

| Scenario | Recommended |
|----------|-------------|
| Phosphorylation using CHARMM36m SP1/SP2 | **Example 05** (fast, no ITP needed) |
| Phosphorylation with custom ITP parameters | **Example 06** (this example) |
| Any other PTM or non-standard AA | **Example 06** |

---

## Detailed walkthrough

## Step 1 — Prepare the wildtype topology

Run `pdb2gmx` with `-ignh` to strip and consistently rebuild hydrogens. The output `.gro` normalises atom names; `wt.top` is discarded.

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

## Step 2 — Prepare the NSAA force-field patch

`prepare_nsaa_ff` reads the methyl-capped ITP fragment, builds a proper RTP entry for SEP by combining the Ser backbone with SEP sidechain parameters, and writes the cross-boundary bonded supplement needed by `pmx gentop`.

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

---

## Step 3 — Build the hybrid structure

`pmx mutate_nsaa` replaces Ser42 with the `S2SEP` hybrid (state A: Ser sidechain; state B: SEP phosphate group) and writes `S2SEP.mtp` for `pmx gentop`. Adjust `--resid 42` to the pmx-renumbered position of your target Ser.

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

---

## Step 4 — Build the hybrid topology

Do **not** pass `-ignh` — `pmx mutate_nsaa` has already positioned all atoms. `pdb2gmx` reads the `SEP` residue definition from `SEP.rtp` placed in the FF directory during Step 2.

```bash
gmx pdb2gmx \
    -f      mutant.pdb \
    -o      processed.gro \
    -p      topol.top \
    -ff     charmm36m-mut \
    -water  tip3p
```

---

## Step 5 — Fill B-state bonded terms

`pmx gentop` uses `S2SEP.mtp` to locate the hybrid residue and `SEP_bonded.itp` to supply the GAFF/cross-boundary bonded parameters for the phosphate group in state B. State B carries two extra negative charges.

```bash
pmx gentop \
    -p                  topol.top \
    -o                  pmxtop.top \
    -ff                 charmm36m-mut \
    --extra_mtp         S2SEP.mtp \
    --supplement_bonded SEP_bonded.itp
```

Expected output:
```
log_> Hybrid Residue -> 42 | S2SEP
log_> Making bonds for state B -> ...
log_> Total charge of state A = 0
log_> Total charge of state B = -2
```

---

## Step 6 — Solvate and add ions

```bash
gmx solvate -cp processed.gro -cs spc216.gro -p pmxtop.top -o solvated.gro

gmx grompp -f mdp/em.mdp -c solvated.gro -r solvated.gro \
           -p pmxtop.top -o ions.tpr -maxwarn 1
printf '13\n' | gmx genion -s ions.tpr -pname K -nname CL \
           -neutral -o ions.gro -p pmxtop.top
```

---

## Step 7 — Energy minimise

```bash
gmx grompp -f mdp/em.mdp -c ions.gro -r ions.gro \
           -p pmxtop.top -o em.tpr -maxwarn 1
gmx mdrun -v -deffnm em -ntmpi 1
```

---

## Step 8 — Equilibrate and run endpoint simulations

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

`results.txt` gives ΔG1 (Ser→SEP free energy cost in the protein). Repeat Steps 1–10 with a short Ser-containing reference peptide in water to get ΔG2, then ΔΔG = ΔG1 − ΔG2.

```bash
pmx analyze \
    -fA transitionA/frame*/ti.xvg \
    -fB transitionB/frame*/ti.xvg
```

```
ΔΔG = ΔG1 − ΔG2
```
