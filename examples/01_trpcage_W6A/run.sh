#!/usr/bin/env bash
# =============================================================================
# Example 1 — Trp Cage W6A: standard pmx amino acid substitution
#
# System  : 1L2Y Trp Cage miniprotein
# Mutation: Trp6 -> Ala  (hybrid residue W2A)
# FF      : charmm36m-mut
# =============================================================================
set -euo pipefail

FF=charmm36m-mut
WATER=tip3p
INPUT=input/1L2Y.pdb

# ── 0. Fetch input if not present ──────────────────────────────────────────
if [[ ! -f "$INPUT" ]]; then
    echo ">>> Fetching 1L2Y from RCSB ..."
    bash input/fetch_input.sh
fi

# ── 1. Set the GMXLIB path so pmx finds its force-field data ───────────────
eval "$(python3 -c 'from pmx.gmx import set_gmxlib; import os; set_gmxlib(); print("export GMXLIB="+os.environ["GMXLIB"])')"

# ── 2. Build hybrid structure ───────────────────────────────────────────────
# --resid 6      : Trp is the 6th residue (pmx renumbers from 1)
# --resname A    : one-letter code of the TARGET residue (Ala)
# -ff            : force field name
echo ">>> pmx mutate ..."
pmx mutate \
    -f  "$INPUT" \
    -o  mutant.pdb \
    --resid   6 \
    --resname A \
    -ff "$FF"

# ── 3. Generate standard GROMACS topology ──────────────────────────────────
# pdb2gmx will recognise the W2A hybrid residue name from mutres.rtp
echo ">>> gmx pdb2gmx ..."
gmx pdb2gmx \
    -f      mutant.pdb \
    -o      processed.gro \
    -p      topol.top \
    -ff     "$FF" \
    -water  "$WATER" \
    -ignh

# ── 4. Fill hybrid B states ─────────────────────────────────────────────────
echo ">>> pmx gentop ..."
pmx gentop \
    -p topol.top \
    -o pmxtop.top

echo ""
echo "=== Done ==="
echo "Hybrid structure : mutant.pdb"
echo "Hybrid topology  : pmxtop.top"
echo ""
echo "Next steps: solvate, add ions, energy minimise, equilibrate,"
echo "then run stateA and stateB production simulations."
