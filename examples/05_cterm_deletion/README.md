# Alchemical C-terminal residue deletion: HSP90 MEEVD peptide (1ELR)

```
State A — full peptide                 State B — C-terminally truncated

  M(1)–E(2)–E(3)–V(4)–D(5)               M(1)–E(2)–E(3)–V(4)–[D(5)]
         |                                              |
        ΔG1 (protein)                                ΔG1
         |                                              |
  M(1)–E(2)–E(3)–V(4)–D(5) ─── ΔG2 ─── M(1)–E(2)–E(3)–V(4)–[D(5)]
          (reference peptide in water)

  ΔΔG = ΔG1 − ΔG2      [D] = non-interacting dummy in state B
```

We compute the free energy cost of removing the C-terminal aspartate from
the HSP90 C-terminal MEEVD peptide (PDB: 1ELR chain B).  The same
calculation in a reference peptide in water gives ΔG2; the difference ΔΔG
reflects how the protein environment modulates the cost of the truncation.

**Terminal deletion** is a special pmx mutation type:
- Asp5 stays physically present but becomes fully non-interacting in state B
  (`DUM_` atom types, zero charge and LJ ε)
- The neighbouring residue Val4 gains a dummy oxygen (`DOT2`, type `DUM_OC`,
  charge 0) that becomes the real C-terminal `OT2` carboxylate oxygen in
  state B

The hybrid residue placed on Val4 is called **`VdeC`** (VAL as C-terminal
deletion neighbour).

> **Note on charge states:** Asp5 at the C-terminus bears a formal charge of
> −2 (sidechain COO⁻ plus the C-terminal carboxylate patch).  When it is
> dummified, both charges are zeroed, while VdeC in state B only gains one
> extra carboxylate oxygen (charge −1).  As a result state A = −3 e and
> state B = −2 e.  This is correct physics — the charge asymmetry reflects
> the removal of a charged residue — and must be compensated by a
> counter-ion correction term in the analysis (add/remove one K⁺/Cl⁻ ion
> between the two boxes, or apply a Poisson–Boltzmann charge correction).

---

## Prerequisites — generate terminal-deletion parameters

The `XdeN` and `XdeC` hybrid entries are not part of the standard mutres
files; they are generated once per force field:

```bash
eval "$(python3 -c 'from pmx.gmx import set_gmxlib; import os; set_gmxlib(); print("export GMXLIB="+os.environ["GMXLIB"])')"
FFDIR=$(python3 -c "import os; from pmx.gmx import set_gmxlib; set_gmxlib(); \
    print(os.environ['GMXLIB']+'/charmm36m-mut.ff')")

python3 -m pmx.scripts.generate_term_deletion \
    --ffdir "$FFDIR" \
    --outdir .
```

This writes `mutres_term.rtp` and `mutres_term.mtp` into the current
directory.  `run.sh` copies them into the FF directory automatically, but
you can also write directly to `$FFDIR` if you have write access.

---

## Step 1 — Prepare the wildtype topology

```bash
eval "$(python3 -c 'from pmx.gmx import set_gmxlib; import os; set_gmxlib(); print("export GMXLIB="+os.environ["GMXLIB"])')"

gmx pdb2gmx \
    -f      input/1ELR_peptide.pdb \
    -o      wt.gro \
    -p      wt.top \
    -ff     charmm36m-mut \
    -water  tip3p \
    -ignh
```

This normalises atom names (especially the backbone amide hydrogen) using
GROMACS's internal rename tables.  The output `.gro` is passed to `pmx mutate`.

---

## Step 2 — Build the hybrid structure

Target the residue to be **deleted** with the code `DEL`.  `pmx` detects that
Asp5 is the C-terminus and automatically:
- renames Val4 → `VdeC`
- adds a dummy oxygen (`DOT2`) to Val4 representing the future C-terminal
  `OT2` carboxylate oxygen

```bash
printf "5 DEL\n" | pmx mutate \
    -f      wt.gro \
    -o      mutant.pdb \
    -ff     charmm36m-mut
```

Asp5 itself is left physically in the structure; it will be dummified at the
topology level by `pmx gentop`.

---

## Step 3 — Build the hybrid topology

```bash
# Copy hybrid RTP/MTP into the FF directory so pdb2gmx can find VdeC
FFDIR=$(python3 -c "import os; from pmx.gmx import set_gmxlib; set_gmxlib(); \
    print(os.environ['GMXLIB']+'/charmm36m-mut.ff')")
cp mutres_term.rtp "$FFDIR/"
cp mutres_term.mtp "$FFDIR/"

gmx pdb2gmx \
    -f      mutant.pdb \
    -o      processed.gro \
    -p      topol.top \
    -ff     charmm36m-mut \
    -water  tip3p
```

Do **not** pass `-ignh` — the dummy atom `DOT2` placed by `pmx mutate` must
be preserved; GROMACS cannot rebuild it from the hydrogen database.

`pdb2gmx` will find the `VdeC` residue definition in `mutres_term.rtp` and
generate a standard GROMACS topology that has Asp5 as a regular ASP (to be
dummified later) and Val4 as a `VdeC` interior residue.

---

## Step 4 — Fill B-state bonded terms

```bash
pmx gentop \
    -p          topol.top \
    -o          pmxtop.top \
    -ff         charmm36m-mut \
    --extra_mtp mutres_term.mtp
```

`pmx gentop` will:
- Find the `VdeC` hybrid residue
- Dummify Asp5 (set all atom types to `DUM_*` in state B, charge 0)
- Apply the C-terminal patch to Val4 in state B: C type C→CC, O type O→OC
  (renamed OT1), DOT2→OC as OT2
- Zero the B-state force constants on dihedrals that cross the Val4–Asp5
  boundary (since Asp5 is absent in state B)

Expected output:
```
log_> Hybrid Residue -> 4 | VdeC
log_> Dummifying deletion target: 5 | ASP
log_> Making bonds for state B -> 30 bonds with perturbed atoms
log_> Making angles for state B -> 56 angles with perturbed atoms
log_> Making dihedrals for state B -> 49 dihedrals with perturbed atoms
log_> Removed 2 fake dihedrals
log_> Total charge of state A = -3
log_> Total charge of state B = -2
log_> Decoupling cross-residue dihedrals: VdeC (res 4) -- ASP (res 5)
log_> Zeroed B-state for 13 dihedrals
```

The A/B charge difference (−3 vs −2) is expected: see the note above.

---

## Step 5 — Solvate and add ions

```bash
gmx solvate -cp processed.gro -cs spc216.gro -p pmxtop.top -o solvated.gro

gmx grompp -f mdp/em.mdp -c solvated.gro -r solvated.gro \
           -p pmxtop.top -o ions.tpr -maxwarn 1
printf '13\n' | gmx genion -s ions.tpr -pname K -nname CL \
           -neutral -o ions.gro -p pmxtop.top
```

Because state A carries one more negative charge than state B, you must run
separate ion sets for each endpoint (or apply a charge-correction in
post-processing).  A common approach is to use the same ion configuration and
apply a Born/Poisson–Boltzmann correction at analysis time.

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

`results.txt` gives ΔG1 (free energy of C-terminal Asp deletion in the
protein context).  Repeat with an Asp–Val reference dipeptide (or the
full MEEVD peptide in water) to get ΔG2, then:

```
ΔΔG = ΔG1 − ΔG2
```

A negative ΔΔG indicates Asp5 is thermodynamically less favourable in the
protein than in water (i.e. C-terminal truncation is favoured by the protein
environment).

---

## Notes

- The `VdeC` hybrid residue type is specific to VAL as the C-terminal
  deletion neighbour.  For other residues, the analogous types are `AdeC`
  (ALA), `EdeC` (GLU), etc. — all generated by `generate_term_deletion.py`.
- The backbone amide hydrogen in `XdeN`/`XdeC` hybrid residues is named `HN`
  (standard CHARMM36 convention).  The paired file `mutres_term.arn` contains
  the rule `* H HN`, which instructs pdb2gmx to rename `H` back to `HN` for
  all residues defined in `mutres_term.rtp`.
- Deleting a **charged** C-terminal residue causes a charge asymmetry between
  state A and state B.  Apply a charge-correction term in post-processing
  (e.g. analytical Born correction or separate ion run).
- The complementary example **04_nterm_deletion** demonstrates deleting the
  N-terminal MET1 residue from the same MEEVD peptide.
