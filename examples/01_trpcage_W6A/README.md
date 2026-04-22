# Alchemical amino acid substitution: Trp6→Ala in Trp Cage (1L2Y)

```
State A — wildtype                     State B — mutant

  Trp6 (folded Trp Cage)                 Ala6 (folded Trp Cage)
       |                                        |
      ΔG1 (protein)                           ΔG1
       |                                        |
  Trp6 (unfolded / water)  ─── ΔG2 ───  Ala6 (unfolded / water)

  ΔΔG = ΔG1 − ΔG2
```

We compute the free energy difference between the wildtype (Trp6) and mutant
(Ala6) forms of the Trp Cage miniprotein (PDB: 1L2Y).  Repeating the
calculation in water (unfolded reference leg) gives ΔG2; the difference ΔΔG
is the change in folding stability due to the substitution.

The Trp residue at position 6 is the dominant hydrophobic anchor of the Trp
Cage hydrophobic core; the W6A mutation is experimentally destabilising by
~2 kcal/mol and serves as the textbook pmx demonstration case.

We use the **`W2A`** hybrid residue: state A carries the full Trp side chain,
state B carries a dummy Trp side chain plus the (much shorter) Ala `CB`
position.

---

## Step 1 — Prepare the wildtype topology

```bash
eval "$(python3 -c 'from pmx.gmx import set_gmxlib; import os; set_gmxlib(); print("export GMXLIB="+os.environ["GMXLIB"])')"

gmx pdb2gmx \
    -f      input/1L2Y.pdb \
    -o      wt.gro \
    -p      wt.top \
    -ff     charmm36m-mut \
    -water  tip3p \
    -ignh
```

`-ignh` strips all existing hydrogens so GROMACS rebuilds them consistently
from the hydrogen database.  The wildtype `.gro` is used only to normalise
atom names before `pmx mutate`; the topology (`wt.top`) is discarded.

---

## Step 2 — Build the hybrid structure

```bash
printf "6 A\n" | pmx mutate \
    -f      wt.gro \
    -o      mutant.pdb \
    -ff     charmm36m-mut
```

`pmx mutate` replaces Trp6 with the `W2A` hybrid: the full Trp side chain is
kept as state A atoms while a set of dummy Ala atoms is added for state B.
The mutation code `A` is the one-letter code of the target residue.

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
and dummy atoms.  `pdb2gmx` reads the `W2A` entry from `mutres.rtp` and
generates a standard topology that includes all hybrid atoms.

---

## Step 4 — Fill B-state bonded terms

```bash
pmx gentop \
    -p  topol.top \
    -o  pmxtop.top \
    -ff charmm36m-mut
```

`pmx gentop` reads the `W2A` entry from `mutres.mtp` and fills in the B-state
atom types, charges, and bonded parameters.  Dummy atoms get type `DUM_*`,
charge 0, and ε = 0.

Expected output:
```
log_> Hybrid Residue -> 6 | W2A
log_> Making bonds for state B -> ...
log_> Total charge of state A = 0
log_> Total charge of state B = 0
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

Allow at least 10 ns equilibration per endpoint before harvesting transition
frames.  The Trp Cage is small (20 residues) and equilibrates quickly.

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

`results.txt` gives ΔG1 (Trp→Ala free energy cost in the protein).  Repeat
Steps 1–9 with a short Trp-containing reference peptide in water to get ΔG2,
then:

```
ΔΔG = ΔG1 − ΔG2
```

A positive ΔΔG means the mutation destabilises the folded protein relative to
the unfolded/solvated reference.  The experimentally measured value for W6A
Trp Cage is approximately +2 kcal/mol.

---

## Notes

- pmx renumbers residues from 1 by default, so `6` refers to the 6th residue
  of the chain regardless of PDB numbering.
- 1L2Y is an NMR ensemble; only MODEL 1 (chain A) is used.  `fetch_input.sh`
  extracts it automatically.
- Both state A and state B are charge-neutral; no charge correction is needed.
- For production-quality ΔΔG, use at least 100 transition trajectories per
  direction and converge each endpoint simulation.
