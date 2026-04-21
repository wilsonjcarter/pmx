#!/usr/bin/env bash
# =============================================================================
# Example 3 — Thioredoxin disulfide formation: C2CD pair
#
# System  : 1ERT E. coli thioredoxin (reduced form)
# Mutation: Cys32 + Cys35 -> both C2CD  (disulfide formation hybrid)
# FF      : charmm36m-mut
#
# State A: both cysteines reduced (free -SH thiols, SG type S)
# State B: Cys32–Cys35 disulfide bond formed (-S-S-, SG type SM)
#
# The SG-SG cross-residue bond + flanking angles/dihedrals are added
# automatically by pmx gentop via _add_disulfide_bonded_terms().
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

# ── 2a. Mutate Cys32 -> C2CD ───────────────────────────────────────────────
# The target one-letter code for the disulfide hybrid state is 'CD'
# (encodes the C2CD hybrid; 'C' = source CYS, 'CD' = CYS_Disulfide target)
echo ">>> pmx mutate: Cys32 -> C2CD ..."
pmx mutate \
    -f  "$INPUT" \
    -o  mutant_step1.pdb \
    --resid   32 \
    --resname CD \
    -ff "$FF"

# ── 2b. Mutate Cys35 -> C2CD ───────────────────────────────────────────────
echo ">>> pmx mutate: Cys35 -> C2CD ..."
pmx mutate \
    -f  mutant_step1.pdb \
    -o  mutant.pdb \
    --resid   35 \
    --resname CD \
    -ff "$FF"

rm -f mutant_step1.pdb

# ── 3. Generate standard GROMACS topology ──────────────────────────────────
# After both pmx mutate steps, Cys32 and Cys35 are named C2CD — they are read
# from mutres.rtp, not as CYS.  pdb2gmx will not see any CYS pair and will not
# attempt automatic disulfide detection.
echo ">>> gmx pdb2gmx ..."
gmx pdb2gmx \
    -f      mutant.pdb \
    -o      processed.gro \
    -p      topol.top \
    -ff     "$FF" \
    -water  "$WATER" \
    -ignh

# ── 4. Fill hybrid B states ─────────────────────────────────────────────────
# gentop:
#   - detects the two C2CD residues (both satisfy is_hybrid())
#   - calls _add_disulfide_bonded_terms() to inject the cross-residue
#     SG32-SG35 bond, CB-SG-SG' and SG-SG'-CB' angles, and proper dihedrals
#   - dummifies HG1 (-> DUM_HS) in state B for both cysteines
echo ">>> pmx gentop ..."
pmx gentop \
    -p topol.top \
    -o pmxtop.top

echo ""
echo "=== Done ==="
echo "Hybrid structure : mutant.pdb    (Cys32 and Cys35 are C2CD)"
echo "Hybrid topology  : pmxtop.top    (SG-SG cross-residue terms included)"
echo ""
echo "State A: reduced thioredoxin  (both Cys free thiols)"
echo "State B: oxidised thioredoxin (Cys32-Cys35 disulfide)"
echo ""
echo "Tip: the SG32-SG35 distance in 1ERT is ~3.9 A (reduced form)."
echo "     Allow generous equilibration before production."
echo "     Run the same transformation with a Cys-Cys dipeptide in water"
echo "     for the reference leg of the thermodynamic cycle."
