#!/usr/bin/env bash
# =============================================================================
# Example 3 — Thioredoxin disulfide formation: Cys32-Cys35
#
# System  : 1ERT E. coli thioredoxin (reduced form)
# Mutation: Cys32 + Cys35 -> CYS2  (disulfide formation)
# FF      : charmm36m-mut
#
# State A: both cysteines reduced (free -SH thiols)
# State B: Cys32-Cys35 disulfide bond formed (-S-S-)
#
# Bond formation is a GENERATIVE operation, not a pre-built library hybrid:
#   pmx mutate --react SG   splits each cysteine so SG is clawed back from the
#                           dummified side chain and survives into state B
#   pmx gentop  --bond      morphs the cross-residue SG-SG bond and its
#                           flanking angles/dihedrals from k=0 (A) to full (B)
# This replaces the retired pre-built C2CD hybrid.
# =============================================================================
set -euo pipefail

FF=charmm36m-mut
WATER=tip3p
INPUT=input/1ERT.pdb

# ── 0. Fetch input if not present ──────────────────────────────────────────
if [[ ! -f "$INPUT" ]]; then
    echo ">>> Fetching 1ERT from RCSB ..."
    bash input/fetch_input.sh
fi

# ── 1. Set GMXLIB ──────────────────────────────────────────────────────────
eval "$(python3 -c 'from pmx.gmx import set_gmxlib; import os; set_gmxlib(); print("export GMXLIB="+os.environ["GMXLIB"])')"

# ── 2. Normalise atom names with pdb2gmx ───────────────────────────────────
echo ">>> gmx pdb2gmx (wildtype — atom-name normalisation) ..."
gmx pdb2gmx \
    -f      "$INPUT" \
    -o      wt.gro \
    -p      wt.top \
    -ff     "$FF" \
    -water  "$WATER" \
    -ignh

# ── 3. Split both cysteines, clawing back SG ───────────────────────────────
# CYS2 is the bound-state (cystine) residue shipped in merged.rtp.  --react names
# the atom that must survive the side-chain decoupling, because it carries the
# bond being formed; both sites are done in ONE call so the two hybrids are
# built against each other.
echo ">>> pmx mutate: Cys32 + Cys35 -> CYS2 (reacting atom SG) ..."
printf "32 CYS2\n35 CYS2\n" > mut.txt
pmx mutate \
    -f      wt.gro \
    -o      mutant.pdb \
    --script mut.txt \
    --react  SG,SG \
    -ff     "$FF"
rm -f mut.txt wt.gro wt.top

# ── 4. Generate GROMACS topology for the hybrid structure ──────────────────
# After both pmx mutate steps, Cys32 and Cys35 are named C2CD — they are read
# from mutres.rtp, not as CYS.  pdb2gmx will not see any CYS pair and will not
# attempt automatic disulfide detection.
# Do NOT pass -ignh: mutant.pdb already has all H atoms from pmx mutate.
echo ">>> gmx pdb2gmx (mutant) ..."
gmx pdb2gmx \
    -f      mutant.pdb \
    -o      processed.gro \
    -p      topol.top \
    -ff     "$FF" \
    -water  "$WATER"

# ── 5. Fill hybrid B states ─────────────────────────────────────────────────
# gentop:
#   - fills the B states of both split cysteines
#   - --bond adds the cross-residue SG32-SG35 bond and its flanking angles and
#     dihedrals as a morph: k=0 in state A (no bond), full value in state B
#   - dummifies the thiol protons in state B for both cysteines
echo ">>> pmx gentop ..."
pmx gentop \
    -p    topol.top \
    -o    pmxtop.top \
    -ff   "$FF" \
    --bond "32:SG-35:SG"

# ── 6. Define box, solvate, and add ions ───────────────────────────────────
# editconf must be called before solvate to define the box vectors.
echo ">>> gmx editconf ..."
gmx editconf \
    -f  processed.gro \
    -o  boxed.gro \
    -bt dodecahedron \
    -d  1.2

echo ">>> gmx solvate ..."
gmx solvate \
    -cp boxed.gro \
    -cs spc216.gro \
    -p  pmxtop.top \
    -o  solvated.gro

echo ">>> gmx genion (0.15 M KCl, neutral) ..."
gmx grompp -f mdp/em.mdp -c solvated.gro -r solvated.gro \
           -p pmxtop.top -o ions.tpr -maxwarn 1
printf '13\n' | gmx genion \
    -s    ions.tpr \
    -pname K  -nname CL \
    -neutral  -conc 0.15 \
    -o    ions.gro \
    -p    pmxtop.top

echo ""
echo "=== Done ==="
echo "Hybrid structure : mutant.pdb    (Cys32 and Cys35 are C2CD)"
echo "Hybrid topology  : pmxtop.top    (SG-SG cross-residue terms included)"
echo "Solvated system  : ions.gro      (ready for energy minimisation)"
echo ""
echo "State A: reduced thioredoxin  (both Cys free thiols)"
echo "State B: oxidised thioredoxin (Cys32-Cys35 disulfide)"
echo ""
echo "Tip: the SG32-SG35 distance in 1ERT is ~3.9 A (reduced form)."
echo "     Allow generous equilibration before production."
echo "     Run the same transformation with a Cys-Cys dipeptide in water"
echo "     for the reference leg of the thermodynamic cycle."
