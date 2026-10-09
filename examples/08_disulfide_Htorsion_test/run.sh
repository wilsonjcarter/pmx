#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Disulfide H-torsion test: does the graph rule's extra HB-CB-SG-SG' torsions
# (which the old disulfide hand-list omitted) actually change the free energy?
#
# Model: dimethyl disulfide (CH3-S-S-CH3), vacuum.  Alchemically turn the S-S
# bond + its induced angles/dihedrals from k=0 (state A) to full GAFF (state B),
# TWO ways:
#   * full      : the complete graph set (induced_bond_morph)  -> incl. H-C-S-S
#   * heavy_only : heavy-atom terms only (reproduces the old hand-list)
# Run reverse non-equilibrium switches for each and compare dG.  The difference
# is the free-energy contribution of the previously-omitted H-torsions.
#
# Requires: antechamber/parmchk2/tleap + parmed, obabel, gmx, pmx (PMX_SRC).
# ---------------------------------------------------------------------------
set -e
SMILES=${SMILES:-"CSSC"}
BOND_A=${BOND_A:-S}; BOND_B=${BOND_B:-S1}
NREP=${NREP:-10}
NSTEPS_EQ=${NSTEPS_EQ:-500000}; NSTOUT=${NSTOUT:-25000}
NSTEPS_SW=${NSTEPS_SW:-100000}
PMX_SRC=${PMX_SRC:-$HOME/Documents/pmx/src}; export PMX_SRC
GMX=${GMX:-gmx}

echo "== 1. dimethyl disulfide, ABCG2/GAFF2 =="
obabel -:"$SMILES" --gen3d -O lig.mol2 2>/dev/null
antechamber -i lig.mol2 -fi mol2 -o lig_abcg2.mol2 -fo mol2 -c abcg2 -nc 0 -at gaff2 -s 0 >ac.log 2>&1
parmchk2 -i lig_abcg2.mol2 -f mol2 -o lig.frcmod -s gaff2
printf 'source leaprc.gaff2\nm = loadmol2 lig_abcg2.mol2\nloadamberparams lig.frcmod\nsaveamberparm m lig.prmtop lig.inpcrd\nquit\n' > leap.in
tleap -f leap.in >leap.log 2>&1
python -c "import parmed as p; a=p.load_file('lig.prmtop',xyz='lig.inpcrd'); a.save('state_b.top',format='gromacs',overwrite=True); a.save('lig.gro',overwrite=True)"

run_variant () {   # $1 = tag (full|heavy_only)
  tag=$1; extra=""; [ "$tag" = "heavy_only" ] && extra="heavy_only"
  echo "== variant: $tag =="
  python build_hybrid.py state_b.top hybrid_$tag.top "$BOND_A" "$BOND_B" $extra
  $GMX editconf -f lig.gro -o box.gro -c -bt cubic -box 4 4 4 >/dev/null 2>&1
  $GMX grompp -f mdp/em.mdp -c box.gro -p hybrid_$tag.top -o em_$tag.tpr -maxwarn 5 >em_$tag.log 2>&1
  $GMX mdrun -deffnm em_$tag -nt 1 >>em_$tag.log 2>&1
  sed "s/NSTEPS_EQ/$NSTEPS_EQ/;s/NSTOUT/$NSTOUT/" mdp/equil.mdp > equil_$tag.mdp
  $GMX grompp -f equil_$tag.mdp -c em_$tag.gro -p hybrid_$tag.top -o equil_$tag.tpr -maxwarn 5 >equil_$tag.log 2>&1
  $GMX mdrun -deffnm equil_$tag -nt 1 >>equil_$tag.log 2>&1
  DL=$(python -c "print(-1.0/$NSTEPS_SW)")
  rm -f works_$tag.txt
  for i in $(seq 1 $NREP); do
    t=$(python -c "print(int($i*$NSTEPS_EQ/($NREP+1)*0.001))")
    echo "System" | $GMX trjconv -s equil_$tag.tpr -f equil_$tag.xtc -o s_${tag}_$i.gro -dump $t >/dev/null 2>&1
    sed "s/NSTEPS_SW/$NSTEPS_SW/;s/DLAMBDA/$DL/;s/SEED/$((1000+i))/" mdp/switch.mdp > sw_${tag}_$i.mdp
    $GMX grompp -f sw_${tag}_$i.mdp -c s_${tag}_$i.gro -p hybrid_$tag.top -o sw_${tag}_$i.tpr -maxwarn 5 >sw_${tag}_$i.log 2>&1
    $GMX mdrun -deffnm sw_${tag}_$i -nt 1 >>sw_${tag}_$i.log 2>&1
    python - "sw_${tag}_$i.xvg" "works_$tag.txt" <<'PY'
import sys, numpy as np
d=np.array([[float(x) for x in l.split()] for l in open(sys.argv[1]) if l[0] not in "@#"])
t,g=d[:,0],d[:,1]; lam=1.0-(t-t[0])/(t[-1]-t[0])
open(sys.argv[2],"a").write(f"{getattr(np,'trapezoid',np.trapz)(g,lam):.4f}\n")
PY
  done
}

run_variant full
run_variant heavy_only

echo "== dG comparison (Jarzynski, reverse) =="
python - <<'PY'
import numpy as np
kT = 8.314462618e-3 * 300.0
# W are reverse works (B->A);  dG(A->B) = +kT ln <exp(-W_rev/kT)>
dGAB = lambda W: kT * np.log(np.mean(np.exp(-W / kT)))
res = {}
for tag in ("full", "heavy_only"):
    W = np.loadtxt(f"works_{tag}.txt", ndmin=1)
    res[tag] = dGAB(W)
    sd = W.std(ddof=1) if len(W) > 1 else 0.0
    print(f"  {tag:10s}: dG(A->B) = {res[tag]:8.2f} kJ/mol   "
          f"(<W_rev> {W.mean():7.2f} +/- {sd:.2f} kJ/mol)")
ddg = res["full"] - res["heavy_only"]
print(f"\n  ==> H-torsion contribution  ddG = {ddg:+.2f} kJ/mol = {ddg/4.184:+.3f} kcal/mol")
print("      (the free-energy effect the old disulfide hand-list omitted)")
PY
