#!/usr/bin/env bash
# =============================================================================
# Example 6 — Serine phosphorylation: Ser → phosphoSer (SP1, monoanionic)
#
# System  : Protein with a target Ser (adjust MUT_RESID below)
# Mutation: Ser -> SP1 (monoanionic O-phosphoserine, charge -1)
# FF      : charmm36m-mut  (SP1 is already parameterised in this FF)
#
# State A: Ser  (-OH,    charge  0)
# State B: SP1  (-OPO3H-, charge -1)
#
# No RDKit, no ITP preparation, no FF patching required.
# For SP2 (dianionic, charge -2): replace P1 with P2 throughout.
#
# Prerequisites:
#   - input/protein.pdb   (protein structure)
#   Run:  bash input/fetch_input.sh   to download the example structure.
# =============================================================================
set -euo pipefail

FF=charmm36m-mut
WATER=tip3p
INPUT=input/protein.pdb

# Set the residue ID to mutate — adjust for your structure!
MUT_RESID=42       # pmx-renumbered position of the target Ser (1-based)
MUT_CODE=P1        # P1 = monoanionic SP1; use P2 for dianionic SP2

# ── 0. Fetch input if not present ──────────────────────────────────────────
if [[ ! -f "$INPUT" ]]; then
    echo ">>> Fetching structure from RCSB ..."
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

# ── 3. Build hybrid structure ───────────────────────────────────────────────
# Mutation code P1 selects the S2P1 hybrid (SER -> monoanionic phosphoSer).
# P2 would select S2P2 (SER -> dianionic phosphoSer).
echo ">>> pmx mutate: Ser${MUT_RESID} -> SP1 ..."
printf "%s %s\n" "$MUT_RESID" "$MUT_CODE" | pmx mutate \
    -f      wt.gro \
    -o      mutant.pdb \
    -ff     "$FF"
rm -f wt.gro wt.top

# ── 4. Generate GROMACS topology for the hybrid structure ──────────────────
# Do NOT pass -ignh: the dummy phosphate atoms placed by pmx mutate must be
# preserved; pdb2gmx would not know how to re-add them from HDB.
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
echo "State A: Ser${MUT_RESID}  (neutral,      charge  0)"
echo "State B: SP1${MUT_RESID}  (monoanionic,  charge -1)"
echo ""
echo "Note: add 1 K+ counter-ion for the state B simulation leg."
