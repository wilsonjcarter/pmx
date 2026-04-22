# Alchemical C-terminal residue deletion: HSP90 MEEVD peptide (1ELR)

<p align="center">
  <img src="../imgs/schematic_p4.png" alt="deletion thermocycle" width="800"/>
</p>

## TLDR

This example computes the free energy cost of truncating a C-terminal residue (Asp5) from the HSP90 MEEVD peptide, making it the only pmx workflow that uses the `DEL` mutation code and the specially generated `XdeC`/`XdeN` hybrid types. The neighbouring Val4 is renamed to `VdeC` and grows a dummy C-terminal oxygen (`DOT2`) that becomes real in state B, while Asp5 is fully dummified — a physically distinct treatment from a standard residue swap.

| | |
|---|---|
| **Hybrid** | `VdeC` (on Val4) |
| **State A** | Val4 + Asp5 (fully interacting); charge −3 |
| **State B** | Val4 as new C-terminus; Asp5 dummified; charge −2 |
| **Charge shift** | −3 → −2 (removing a −2 charged C-terminal residue) |
| **Special steps** | Run `generate_term_deletion` first; copy `mutres_term.rtp/.mtp` to FF dir |

```bash
# Unique mutation step — see run.sh for the full pipeline
printf "5 DEL\n" | pmx mutate \
    -f  wt.gro \
    -o  mutant.pdb \
    -ff charmm36m-mut
```

---

## Detailed walkthrough

## Step 1 — Prepare the wildtype topology

Run `pdb2gmx` to normalise atom names (especially the backbone amide hydrogen) before passing the structure to `pmx mutate`. The output `.gro` is the only product needed from this step; `wt.top` is discarded.

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

---

## Step 2 — Generate terminal-deletion parameters

The `XdeN`/`XdeC` hybrid entries are not part of the standard mutres files and must be generated once per force field. `run.sh` copies the output into the FF directory automatically.

```bash
eval "$(python3 -c 'from pmx.gmx import set_gmxlib; import os; set_gmxlib(); print("export GMXLIB="+os.environ["GMXLIB"])')"
FFDIR=$(python3 -c "import os; from pmx.gmx import set_gmxlib; set_gmxlib(); \
    print(os.environ['GMXLIB']+'/charmm36m-mut.ff')")

python3 -m pmx.scripts.generate_term_deletion \
    --ffdir "$FFDIR" \
    --outdir .
```

---

## Step 3 — Build the hybrid structure

Target the residue to be **deleted** with the code `DEL`. `pmx` detects that Asp5 is the C-terminus and automatically renames Val4 → `VdeC` and adds the dummy oxygen `DOT2`; Asp5 stays physically in the structure and is dummified later.

```bash
printf "5 DEL\n" | pmx mutate \
    -f      wt.gro \
    -o      mutant.pdb \
    -ff     charmm36m-mut
```

---

## Step 4 — Build the hybrid topology

Copy the generated hybrid parameters into the FF directory so `pdb2gmx` can resolve `VdeC`, then process the mutant. Do **not** pass `-ignh` — the dummy atom `DOT2` placed by `pmx mutate` must be preserved.

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

---

## Step 5 — Fill B-state bonded terms

`pmx gentop` locates the `VdeC` hybrid, dummifies Asp5 in state B, applies the C-terminal patch to Val4, and zeroes B-state force constants on dihedrals crossing the Val4–Asp5 boundary. The A/B charge difference (−3 vs −2) is expected and must be accounted for in analysis.

```bash
pmx gentop \
    -p          topol.top \
    -o          pmxtop.top \
    -ff         charmm36m-mut \
    --extra_mtp mutres_term.mtp
```

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

`results.txt` gives ΔG1 (free energy of C-terminal Asp deletion in the protein). Repeat with a reference peptide in water to get ΔG2, then ΔΔG = ΔG1 − ΔG2.

```bash
pmx analyze \
    -fA transitionA/frame*/ti.xvg \
    -fB transitionB/frame*/ti.xvg
```

```
ΔΔG = ΔG1 − ΔG2
```
