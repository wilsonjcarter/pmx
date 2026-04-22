# Alchemical pKa shift: Cys32 CYS→CYM in thioredoxin (1ERT)

<p align="center">
  <img src="../imgs/schematic_p2.png" alt="pka thermocycle" width="800"/>
</p>

We compute the free energy cost of deprotonating Cys32 in *E. coli*
thioredoxin (PDB: 1ERT) and subtract the same quantity measured for a
reference cysteine in water.  The difference ΔΔG gives the shift in the
apparent pKa relative to the intrinsic cysteine pKa (~8.3 in solution).

Thioredoxin's active-site Cys32 has an anomalously low pKa (~6.3) driven
by the local electrostatic environment.  This example quantifies that
shift from first principles.

We use the **`C2CM`** hybrid residue:
- **State A** — protonated: `SG` type `S` (charge −0.23), `HG1` present
  (type `HS`, charge 0.16)
- **State B** — deprotonated: `SG` type `SM` (charge −0.75), `HG1` becomes
  a non-interacting dummy (`DUM_HS`, charge 0)

---

## Step 1 — Prepare the wildtype topology

```bash
eval "$(python3 -c 'from pmx.gmx import set_gmxlib; import os; set_gmxlib(); print("export GMXLIB="+os.environ["GMXLIB"])')"

gmx pdb2gmx \
    -f      input/1ERT.pdb \
    -o      wt.gro \
    -p      wt.top \
    -ff     charmm36m-mut \
    -water  tip3p \
    -ignh
```

`-ignh` strips all existing hydrogens so GROMACS rebuilds them consistently.
The wildtype `.gro` is used only to normalise atom names before `pmx mutate`;
the topology (`wt.top`) is discarded.

---

## Step 2 — Build the hybrid structure

```bash
printf "32 CM\n" | pmx mutate \
    -f      wt.gro \
    -o      mutant.pdb \
    -ff     charmm36m-mut
```

`pmx mutate` replaces Cys32 with the `C2CM` hybrid.  The mutation code `CM`
is the pmx extended one-letter code for CYM (deprotonated cysteine).
State A retains the full `SH` group; state B carries a dummy `HG1`
(`DUM_HS`) so the topology has the same number of atoms throughout.

---

## Step 3 — Build the hybrid topology

```bash
gmx pdb2gmx \
    -f      mutant.pdb \
    -o      processed.gro \
    -p      topol.top \
    -ff     charmm36m-mut \
    -water  tip3p
```

Do **not** pass `-ignh` — `pmx mutate` has already positioned all hydrogens
and dummy atoms.  Because Cys32 is renamed `C2CM` (not `CYS`), `pdb2gmx`
will not attempt to form a disulfide involving it.  The remaining Cys35 has
no partner, so no disulfide is detected.

---

## Step 4 — Fill B-state bonded terms

```bash
pmx gentop \
    -p  topol.top \
    -o  pmxtop.top \
    -ff charmm36m-mut
```

`pmx gentop` reads the `C2CM` entry from `mutres.mtp` and fills in the
B-state atom types, charges, and bonded parameters for Cys32.

Expected output:
```
log_> Hybrid Residue -> 32 | C2CM
log_> Making bonds for state B -> ...
log_> Total charge of state A = -5
log_> Total charge of state B = -6
```

State B carries one extra negative charge (the thiolate).  Add one K⁺
counter-ion to the state B simulation box (or run with a background
neutralising charge and apply an analytical correction).

---

## Step 5 — Solvate and add ions

```bash
gmx solvate -cp processed.gro -cs spc216.gro -p pmxtop.top -o solvated.gro

gmx grompp -f mdp/em.mdp -c solvated.gro -r solvated.gro \
           -p pmxtop.top -o ions.tpr -maxwarn 1
printf '13\n' | gmx genion -s ions.tpr -pname K -nname CL \
           -neutral -o ions.gro -p pmxtop.top
```

Neutralise relative to the **state A** charge.  For the state B production
simulation you will need one additional K⁺ ion (or apply a charge correction
in post-processing).

---

## Step 6 — Energy minimise

```bash
gmx grompp -f mdp/em.mdp -c ions.gro -r ions.gro \
           -p pmxtop.top -o em.tpr -maxwarn 1
gmx mdrun -v -deffnm em -ntmpi 1
```

---

## Step 7 — Equilibrate and run endpoint simulations

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

Allow at least 10 ns equilibration per endpoint.  The protonation state
affects the local hydrogen-bond network around the active site; give the
system enough time to equilibrate these contacts.

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

```bash
pmx analyze \
    -fA transitionA/frame*/ti.xvg \
    -fB transitionB/frame*/ti.xvg
```

`results.txt` gives ΔG1 (deprotonation free energy in the protein).  Repeat
Steps 1–9 with a short Cys-containing reference peptide in water to get ΔG2,
then:

```
ΔΔG = ΔG1 − ΔG2

pKa = pKa_ref + ΔΔG / (RT ln 10)
```

Using pKa_ref = 8.3 (solution Cys), a negative ΔΔG means deprotonation is
more favourable in the protein → lower pKa (Cys32 in thioredoxin gives
pKa ≈ 6.3, so ΔΔG ≈ −2.7 kcal/mol).

---

## Notes

- The `C2CM` hybrid uses a simple single-atom perturbation: only `SG` (type
  change S→SM, charge change −0.23→−0.75) and `HG1` (HS→DUM_HS) are
  perturbed.  All backbone and other side-chain atoms are unchanged.
- State B carries one more negative charge than state A.  A charge-correction
  term must be included for quantitative pKa predictions in periodic boundary
  conditions (see Rocklin *et al.* 2013 or the Born/Poisson correction).
- 1ERT also contains Cys35 at the active site (normal pKa ~8).  This example
  targets Cys32 only; for Cys35 use residue number 35 in `pmx mutate`.
- The same `C2CM` transformation applied to a CYM-parameterised reference
  cysteine (starting from the deprotonated form) gives the reverse leg.
