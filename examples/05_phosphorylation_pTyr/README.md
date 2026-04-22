# Alchemical tyrosine phosphorylation: Tyr → phosphoTyr (YP1) in Lck SH2 domain (1AOT)

<p align="center">
  <img src="../imgs/schematic_p5.png" alt="phospho thermocycle" width="800"/>
</p>

## TLDR

This example computes the phosphorylation ΔΔG (Tyr → monoanionic phosphotyrosine YP1) using CHARMM36m's built-in YP1 residue in the Lck SH2 domain (PDB 1AOT) — no force-field patching required. The hybrid `Y2P1` morphs the tyrosine phenol into a phosphate by converting `HH` to a non-interacting dummy and promoting five dummy phosphate atoms (DP, DO2, DH2, DO3, DO4) to real ones, capturing the electrostatic change that drives SH2-domain phosphorylation recognition.

| | |
|---|---|
| **Hybrid** | `Y2P1` |
| **State A** | Tyr: `OH` type `OH1`, `HH` present; charge 0 |
| **State B** | YP1: `OH` → `ON2B`, `HH` → `DUM_H`, dummy `DP/DO2/DH2/DO3/DO4` → real; charge −1 |
| **Charge shift** | 0 → −1 — **doublebox required** (see below) |
| **Special steps** | None — YP1 is already in `charmm36m-mut`; use `P2` / `Y2P2` for dianionic (−2) form |

```bash
# Unique mutation step — see run.sh for the full pipeline
printf "63 P1\n" | pmx mutate \
    -f  wt.gro \
    -o  mutant.pdb \
    -ff charmm36m-mut
```

---

## Detailed walkthrough

### Step 1 — Fetch the structure and prepare the wildtype topology

Download 1AOT (Lck SH2 domain), then run `pdb2gmx` to normalise atom names. `wt.top` is discarded; only `wt.gro` feeds into `pmx mutate`.

```bash
bash input/fetch_input.sh   # downloads 1AOT, extracts chain A, prints Tyr positions

eval "$(python3 -c 'from pmx.gmx import set_gmxlib; import os; set_gmxlib(); print("export GMXLIB="+os.environ["GMXLIB"])')"

gmx pdb2gmx \
    -f      input/1AOT_A.pdb \
    -o      wt.gro \
    -p      wt.top \
    -ff     charmm36m-mut \
    -water  tip3p \
    -ignh
```

Check the Tyr list printed by `fetch_input.sh` to confirm the pmx-renumbered position of your target residue (the script uses 1-based sequential numbering matching pmx's convention).

---

### Step 2 — Build the hybrid structure

Mutation code `P1` selects the `Y2P1` hybrid (Tyr → monoanionic phosphotyrosine). Replace `63` with the residue number identified in Step 1.

```bash
printf "63 P1\n" | pmx mutate \
    -f      wt.gro \
    -o      mutant.pdb \
    -ff     charmm36m-mut
```

`pmx mutate` changes `OH` from type `OH1` to the dummy-bearing `ON2B`, converts `HH` to a zero-charge dummy, and adds five dummy phosphate atoms around the phenol oxygen. Do **not** pass `-ignh` in subsequent steps.

---

### Step 3 — Build the hybrid topology

`pdb2gmx` reads the `Y2P1` residue definition from `mutres.rtp` in `charmm36m-mut` and generates a standard GROMACS topology.

```bash
gmx pdb2gmx \
    -f      mutant.pdb \
    -o      processed.gro \
    -p      topol.top \
    -ff     charmm36m-mut \
    -water  tip3p
```

---

### Step 4 — Fill B-state bonded terms

`pmx gentop` reads `Y2P1` from `mutres.mtp` and fills in all B-state atom types, charges, and bonded parameters.

```bash
pmx gentop \
    -p  topol.top \
    -o  pmxtop.top \
    -ff charmm36m-mut
```

Expected output:
```
log_> Hybrid Residue -> 63 | Y2P1
log_> Making bonds for state B -> ...
log_> Total charge of state A =  0
log_> Total charge of state B = -1
```

---

### Step 5 — Charge-neutral setup with pmx doublebox

Phosphorylation shifts the system charge from 0 to −1. In a periodic simulation box this
charge change introduces finite-size artefacts that can bias ΔG by several kJ/mol. The
correct approach is the **single-box double-system** method: place the protein system and a
reference peptide in the *same* box so that one gains charge while the other loses it,
keeping the total box charge constant throughout the alchemical transition.

```
                          alchemical λ: 0 → 1
  protein in box:   Tyr  ──────────────────────►  pTyr     Δq = −1
  reference in box: pTyr ──────────────────────►  Tyr      Δq = +1
  ─────────────────────────────────────────────────────────────────
  net charge change in box:                                    0  ✓
```

ΔΔG is recovered directly from the combined work values: `pmx analyze` sees the total
work of both simultaneous transformations. The reference (free peptide in solution) cancels
ΔG(phosphorylation in water), leaving ΔΔG of SH2 binding — how much more tightly the SH2
domain recruits the phosphorylated peptide vs. the unphosphorylated form.

#### 5a — Prepare the reference leg

The reference is the **Tyr-containing phosphopeptide extracted from 1AOT** — the biological
SH2 ligand in free solution, not a minimal ACE-Tyr-NME tripeptide. Using the actual peptide
captures the free energy of phosphorylation in the same sequence context and avoids
artefacts from end-cap chemistry. The combined ΔΔG directly quantifies the binding
preference of the SH2 domain for the phosphorylated vs. unphosphorylated form of the
ligand (see `input/fetch_input.sh` to extract the peptide from 1AOT).

```bash
# Reference: Tyr-peptide extracted from 1AOT (provide as input/1AOT_peptide.pdb)
gmx pdb2gmx \
    -f      input/1AOT_peptide.pdb \
    -o      ref_wt.gro \
    -p      ref_wt.top \
    -ff     charmm36m-mut \
    -water  tip3p \
    -ignh

printf "1 P1\n" | pmx mutate \
    -f      ref_wt.gro \
    -o      ref_mutant.pdb \
    -ff     charmm36m-mut

gmx pdb2gmx \
    -f      ref_mutant.pdb \
    -o      ref_processed.gro \
    -p      ref_topol.top \
    -ff     charmm36m-mut \
    -water  tip3p

pmx gentop \
    -p  ref_topol.top \
    -o  ref_pmxtop.top \
    -ff charmm36m-mut
```

#### 5b — Combine into one box

```bash
pmx doublebox \
    -f1 processed.gro \
    -f2 ref_processed.gro \
    -o  doublebox.gro \
    -r  2.5 \
    -d  1.5
```

This places the protein and reference peptide in a single rectangular box separated by at
least 2.5 nm, with 1.5 nm to the box wall.

#### 5c — Merge topologies, solvate, and add ions

Combine the two topology files by appending the reference molecule section to the protein
topology. Then solvate the combined box and neutralise with `gmx genion`:

```bash
# Solvate the combined box
gmx solvate -cp doublebox.gro -cs spc216.gro \
            -p pmxtop.top -o solvated.gro

# Add ions — system is charge-neutral by construction so no extra counter-ions needed
gmx grompp -f mdp/em.mdp -c solvated.gro -r solvated.gro \
           -p pmxtop.top -o ions.tpr -maxwarn 1
printf '13\n' | gmx genion -s ions.tpr -pname K -nname CL \
           -neutral -conc 0.15 -o ions.gro -p pmxtop.top
```

---

### Step 6 — Energy minimise

```bash
gmx grompp -f mdp/em.mdp -c ions.gro -r ions.gro \
           -p pmxtop.top -o em.tpr -maxwarn 1
gmx mdrun -v -deffnm em -ntmpi 1
```

---

### Step 7 — Equilibrate and run endpoint simulations

Tyrosine phosphorylation introduces a large electrostatic change; allow at least 10 ns equilibration per endpoint before harvesting transition frames.

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

### Step 8 — Non-equilibrium transition simulations

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

### Step 9 — Analyse

Because both legs (protein and reference peptide) run simultaneously in the same box,
`pmx analyze` operates on the combined work values and directly yields ΔΔG:

```bash
pmx analyze \
    -fA transitionA/frame*/ti.xvg \
    -fB transitionB/frame*/ti.xvg
```

The output ΔG is ΔΔG = ΔG(phosphorylation when bound to SH2) − ΔG(phosphorylation of the
free peptide in water). A negative ΔΔG means the SH2 domain specifically stabilises the
phosphorylated peptide — the thermodynamic signature of phospho-recognition and the direct
binding ΔΔG between the pTyr and Tyr forms of the ligand.
