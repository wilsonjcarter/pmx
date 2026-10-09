#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# BTK-style covalent-inhibitor demo: alchemical acrylamide -> Cys thioether
# bond formation, in vacuum, via the pmx covalent framework.
#
#   state B (bonded thioether)  --reverse switch-->  state A (bond off)
#
# The forming S-C bond + its graph-induced angles/dihedrals (selected by
# pmx.covalent.edit.induced_bond_morph) morph k=0 (A) <-> full GAFF (B).
# We sample the well-defined state B and run N reverse non-equilibrium
# switches, then estimate dG by Jarzynski (reverse).
#
# Requires: antechamber/parmchk2/tleap + parmed (conda env `openff` or an
# AmberTools install), obabel, gmx, and pmx on PYTHONPATH (PMX_SRC).
# ---------------------------------------------------------------------------
set -e

# --- knobs (override via env) ---
SMILES=${SMILES:-"CSCCC(=O)NC"}     # Cys-SG mimic + acrylamide Michael adduct
BOND_A=${BOND_A:-S}                 # forming-bond atom names (see mapping note)
BOND_B=${BOND_B:-C1}
NREP=${NREP:-10}                    # number of switches
NSTEPS_EQ=${NSTEPS_EQ:-500000}      # state-B equilibration steps (dt=1fs -> 500 ps)
NSTOUT=${NSTOUT:-50000}             # save a frame every ... (-> NREP frames)
NSTEPS_SW=${NSTEPS_SW:-100000}      # switch length steps (dt=0.5fs -> 50 ps)
PMX_SRC=${PMX_SRC:-$HOME/Documents/pmx/src}
export PMX_SRC
GMX=${GMX:-gmx}

echo "== 1. ligand structure + ABCG2/GAFF2 params =="
obabel -:"$SMILES" --gen3d -O lig.mol2 2>/dev/null
antechamber -i lig.mol2 -fi mol2 -o lig_abcg2.mol2 -fo mol2 -c abcg2 -nc 0 -at gaff2 -s 0 >ac.log 2>&1
parmchk2 -i lig_abcg2.mol2 -f mol2 -o lig.frcmod -s gaff2
cat > leap.in <<EOF
source leaprc.gaff2
m = loadmol2 lig_abcg2.mol2
loadamberparams lig.frcmod
saveamberparm m lig.prmtop lig.inpcrd
quit
EOF
tleap -f leap.in >leap.log 2>&1
python -c "import parmed as p; a=p.load_file('lig.prmtop',xyz='lig.inpcrd'); a.save('state_b.top',format='gromacs',overwrite=True); a.save('lig.gro',overwrite=True)"

echo "== 2. dual-state hybrid via induced_bond_morph =="
python build_hybrid.py state_b.top hybrid.top "$BOND_A" "$BOND_B"

echo "== 3. box + minimize (state B) =="
$GMX editconf -f lig.gro -o box.gro -c -bt cubic -box 4 4 4 >/dev/null 2>&1
$GMX grompp -f mdp/em.mdp -c box.gro -p hybrid.top -o em.tpr -maxwarn 5 >em.log 2>&1
$GMX mdrun -deffnm em >>em.log 2>&1

echo "== 4. equilibrate state B, harvest $NREP frames =="
sed "s/NSTEPS_EQ/$NSTEPS_EQ/;s/NSTOUT/$NSTOUT/" mdp/equil.mdp > equil.run.mdp
$GMX grompp -f equil.run.mdp -c em.gro -p hybrid.top -o equil.tpr -maxwarn 5 >equil.log 2>&1
$GMX mdrun -deffnm equil >>equil.log 2>&1

echo "== 5. $NREP reverse switches (B -> A) =="
DL=$(python -c "print(-1.0/$NSTEPS_SW)")
rm -f works.txt
# extract frames evenly from the equilibration trajectory
for i in $(seq 1 $NREP); do
  t=$(python -c "print(int($i*$NSTEPS_EQ/$NREP*0.001))")   # ps
  echo "System" | $GMX trjconv -s equil.tpr -f equil.xtc -o start_$i.gro -dump $t >/dev/null 2>&1
  sed "s/NSTEPS_SW/$NSTEPS_SW/;s/DLAMBDA/$DL/;s/SEED/$((1000+i))/" mdp/switch.mdp > sw_$i.mdp
  $GMX grompp -f sw_$i.mdp -c start_$i.gro -p hybrid.top -o sw_$i.tpr -maxwarn 5 >sw_$i.log 2>&1
  $GMX mdrun -deffnm sw_$i >>sw_$i.log 2>&1
  python - "$i" <<'PY'
import sys, numpy as np
i = sys.argv[1]
d = np.array([[float(x) for x in l.split()] for l in open(f"sw_{i}.xvg") if l[0] not in "@#"])
t, g = d[:,0], d[:,1]
lam = 1.0 - (t - t[0])/(t[-1]-t[0])     # lambda 1 -> 0
W = getattr(np,"trapezoid",np.trapz)(g, lam)                     # reverse work (B->A)
open("works.txt","a").write(f"{W:.4f}\n")
print(f"  switch {i}: W_rev(B->A) = {W:8.2f} kJ/mol")
PY
done

echo "== 6. Jarzynski dG =="
python analyze.py works.txt
