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

# ── 2. Normalise atom names with pdb2gmx ───────────────────────────────────
# pmx mutate requires CHARMM-convention atom names.  Run pdb2gmx on the raw
# PDB first to get a standardised GRO, then use that as input to pmx mutate.
echo ">>> gmx pdb2gmx (wildtype — atom-name normalisation) ..."
gmx pdb2gmx \
    -f      "$INPUT" \
    -o      wt.gro \
    -p      wt.top \
    -ff     "$FF" \
    -water  "$WATER" \
    -ignh

# ── 3. Build hybrid structure ───────────────────────────────────────────────
# Cys32 is the active-site cysteine with depressed pKa (~6.3).
# CM is the pmx extended one-letter code for CYM (deprotonated Cys).
echo ">>> pmx mutate ..."
printf "32 CM\n" > mut.txt
pmx mutate \
    -f      wt.gro \
    -o      mutant.pdb \
    --script mut.txt \
    -ff     "$FF"
rm -f mut.txt wt.gro wt.top

# ── 4. Generate GROMACS topology for the hybrid structure ──────────────────
# After pmx mutate, Cys32 is renamed to C2CM — pdb2gmx reads it from mutres.rtp
# and will not form any disulfide involving it.  The only remaining CYS is
# Cys35, which has no partner, so no disulfide detection is triggered.
# Do NOT pass -ignh: mutant.pdb already has all H atoms from pmx mutate.
echo ">>> gmx pdb2gmx (mutant) ..."
gmx pdb2gmx \
    -f      mutant.pdb \
    -o      processed.gro \
    -p      topol.top \
    -ff     "$FF" \
    -water  "$WATER"

# ── 5. Fill hybrid B states ─────────────────────────────────────────────────
echo ">>> pmx gentop ..."
pmx gentop \
    -p  topol.top \
    -o  pmxtop.top \
    -ff "$FF"

echo ""
echo "=== Done ==="
echo "Hybrid structure : mutant.pdb"
echo "Hybrid topology  : pmxtop.top"
echo ""
echo "State A charge: $(python3 -c "
from pmx.forcefield import Topology; t=Topology('pmxtop.top'); print('%.0f' % t.get_qA())
") e"
echo "State B charge: $(python3 -c "
from pmx.forcefield import Topology; t=Topology('pmxtop.top'); print('%.0f' % t.get_qB())
") e"
echo ""
echo "Note: add one K+ counter-ion for the state B (charged) leg."
