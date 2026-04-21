#!/usr/bin/env bash
# =============================================================================
# Example 2 — Thioredoxin CYS32→CYM: cysteine protonation state FEP
#
# System  : 1ERT E. coli thioredoxin (reduced form)
# Mutation: Cys32 -> Cym  (hybrid residue C2CM)
# FF      : charmm36m-mut
#
# State A: CYS (protonated, -SH, neutral)
# State B: CYM (deprotonated, -S-, charge -1)
#
# Run the same transformation in water (reference leg) and in the protein
# (protein leg) to get the thermodynamic cycle ΔΔG = ΔpKa × RT ln(10).
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

# ── 2. Build hybrid structure ───────────────────────────────────────────────
# Cys32 is the active-site cysteine with depressed pKa (~6.3).
# --resname CM : one-letter code for CYM in the pmx scheme
echo ">>> pmx mutate ..."
pmx mutate \
    -f  "$INPUT" \
    -o  mutant.pdb \
    --resid   32 \
    --resname CM \
    -ff "$FF"

# ── 3. Generate standard GROMACS topology ──────────────────────────────────
# After pmx mutate, Cys32 is renamed to C2CM — pdb2gmx reads it from mutres.rtp
# and will not form any disulfide involving it.  The only remaining CYS is
# Cys35, which has no partner, so no disulfide detection is needed.
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
echo "State A charge: $(python3 -c \"
from pmx.forcefield import Topology; t=Topology('pmxtop.top'); print('%.0f' % t.get_qA())
\") e"
echo "State B charge: $(python3 -c \"
from pmx.forcefield import Topology; t=Topology('pmxtop.top'); print('%.0f' % t.get_qB())
\") e"
echo ""
echo "Note: add one K+ counter-ion for the state B (charged) leg."
