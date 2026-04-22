# Alchemical phosphorylation: Ser → phosphoSer (SP1) in thioredoxin

<p align="center">
  <img src="../imgs/schematic_p5.png" alt="phospho thermocycle" width="800"/>
</p>

We compute the free energy difference between the unphosphorylated (Ser) and
monoanionic phosphoserine (SP1) forms of a target serine.  Repeating the
calculation in water gives ΔG2; the difference ΔΔG is the effect of the
protein environment on the phosphorylation affinity.

This example uses the **standard `pmx mutate` workflow** with the **`SP1`
(monoanionic phosphoserine, −1 charge)** residue that is already parameterised
in the CHARMM36m force field.  No RDKit, no ITP preparation, and no
`prepare_nsaa_ff` step are required.

We use the **`S2P1`** hybrid residue:
- **State A** — unphosphorylated: `OG` type `OH1` (charge −0.66), `HG1`
  present (type `H`, charge 0.43); CB charge +0.05
- **State B** — phosphorylated SP1: `OG` becomes bridging `ON2` (charge
  −0.62); `HG1` becomes a non-interacting dummy (`DUM_H`); dummy phosphate
  atoms (`DP`, `DO1P`, `DO2P`, `DOT`, `DHT`) become real P/ON3/ON4/HN4

Net charge shift: 0 → −1.

> For dianionic phosphoserine (SP2, −2 charge) at pH 7, use mutation code
> `P2` instead of `P1` everywhere below.  The S2P2 hybrid is constructed
> analogously and fully parameterised in the same force field.

---

## Step 1 — Prepare the wildtype topology

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

`-ignh` strips all existing hydrogens so GROMACS rebuilds them consistently.
The wildtype `.gro` is used only to normalise atom names before `pmx mutate`;
the topology (`wt.top`) is discarded.

---

## Step 2 — Build the hybrid structure

```bash
printf "42 P1\n" | pmx mutate \
    -f      wt.gro \
    -o      mutant.pdb \
    -ff     charmm36m-mut
```

`pmx mutate` replaces Ser42 with the `S2P1` hybrid.  The mutation code `P1`
is the pmx extended one-letter code for SP1 (monoanionic phosphoserine).
You can also supply the three-letter code directly: `printf "42 SP1\n"`.

`pmx mutate` inserts:
- the `OG` atom with its type changed (OH1 → ON2 in state B)
- a dummy `HG1` that becomes non-interacting in state B
- four dummy phosphate heavy atoms (`DP`, `DO1P`, `DO2P`, `DOT`) and one
  dummy hydrogen (`DHT`) that become real in state B

Adjust residue number 42 to the pmx-renumbered position of your target Ser
(pmx renumbers from 1 within each chain by default).

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

Do **not** pass `-ignh` — `pmx mutate` has already positioned all atoms.
`pdb2gmx` reads the `S2P1` residue definition from `mutres.rtp` (part of
`charmm36m-mut`) and generates the standard GROMACS topology.

---

## Step 4 — Fill B-state bonded terms

```bash
pmx gentop \
    -p  topol.top \
    -o  pmxtop.top \
    -ff charmm36m-mut
```

`pmx gentop` reads the `S2P1` entry from `mutres.mtp` and fills in the
B-state atom types, charges, and bonded parameters.

Expected output:
```
log_> Hybrid Residue -> 42 | S2P1
log_> Making bonds for state B -> ...
log_> Total charge of state A = 0
log_> Total charge of state B = -1
```

State B carries one extra negative charge.  One additional K⁺ counter-ion
is required for the state B simulation leg.

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

Phosphorylation introduces large electrostatic changes; allow at least 10 ns
equilibration per endpoint.  Monitor the phosphate–protein contacts and the
local hydrogen-bond network to confirm the environment has stabilised before
harvesting transition frames.

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

`results.txt` gives ΔG1 (Ser→SP1 free energy cost in the protein).  Repeat
Steps 1–9 with a short Ser-containing reference peptide in water to get ΔG2,
then:

```
ΔΔG = ΔG1 − ΔG2
```

A negative ΔΔG means phosphorylation is more favourable in the protein
environment than in water (i.e. the protein stabilises the phosphorylated
form).

---

## Notes

- The `S2P1` hybrid uses the CHARMM36m `SP1` residue that is already present
  in `merged.rtp`.  No force-field patching or ITP preparation is required.
  For a fully de-novo non-standard amino acid workflow (custom residues not
  in the FF library) see `examples/07_nsaa_denovo/`.
- Mutation code `P2` / three-letter code `SP2` selects the **dianionic**
  (−2) form.  The state B charge becomes −2; add two K⁺ counter-ions for
  state B instead of one.
- The charge difference between state A (0) and state B (−1) requires a
  charge-correction term for quantitative free energies in periodic boundary
  conditions.  See Rocklin *et al.* (2013) or the Born/Poisson correction.
- pmx renumbers residues from 1 within each chain by default.  Adjust the
  residue number in Step 2 to match your target serine.
