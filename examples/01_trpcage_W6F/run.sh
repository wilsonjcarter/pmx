#!/usr/bin/env bash
# =============================================================================
# Example 1 — Trp Cage W6F: standard pmx amino acid substitution
#
# System  : 1L2Y Trp Cage miniprotein
# Mutation: Trp→Phe  (hybrid residue W2F)
# FF      : charmm36m-mut
#
# Workflow
# --------
# 1. gmx pdb2gmx on the raw PDB  ->  wt.gro  (standardises atom names)
# 2. pmx mutate   on wt.gro      ->  mutant.pdb  (inserts W2F hybrid)
# 3. gmx pdb2gmx  on mutant.pdb  ->  topol.top  (builds topology for W2F)
# 4. pmx gentop   on topol.top   ->  pmxtop.top  (fills B states)
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

# ── 2. Normalise atom names with pdb2gmx ───────────────────────────────────
# pmx mutate expects CHARMM-convention names (HN, HB1, HB2 …) rather than
# the names typically found in raw PDB files (H, HB2, HB3 …).
# Running pdb2gmx first on the wildtype PDB produces a properly named GRO.
echo ">>> gmx pdb2gmx (wildtype — atom-name normalisation) ..."
gmx pdb2gmx \
    -f      "$INPUT" \
    -o      wt.gro \
    -p      wt.top \
    -ff     "$FF" \
    -water  "$WATER" \
    -ignh

# ── 3. Build hybrid structure ───────────────────────────────────────────────
# pmx mutate uses a --script file: "residue_id target_one_letter" per line.
# Residue 6 is Trp→Phe).
echo ">>> pmx mutate ..."
printf "6 F\n" > mut.txt
pmx mutate \
    -f      wt.gro \
    -o      mutant.pdb \
    --script mut.txt \
    -ff     "$FF"
rm -f mut.txt wt.gro wt.top

# ── 4. Generate GROMACS topology for the hybrid structure ──────────────────
# mutant.pdb already has all H atoms placed by pmx mutate; do NOT pass -ignh.
# pdb2gmx recognises the W2F hybrid residue from mutres.rtp.
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
echo "Next steps: solvate, add ions, energy minimise, equilibrate,"
echo "then run stateA and stateB production simulations."
