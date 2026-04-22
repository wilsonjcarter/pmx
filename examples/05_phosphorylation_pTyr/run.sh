#!/usr/bin/env bash
# =============================================================================
# Example 5 — Tyrosine phosphorylation: Tyr → phosphoTyr (TP1, monoanionic)
#
# System  : Lck SH2 domain (PDB 1AOT, chain A)
# Mutation: Tyr -> TP1 (monoanionic O-phosphotyrosine, charge -1)
# FF      : charmm36m-mut  (TP1 is already parameterised in this FF)
#
# State A: Tyr  (-OH,      charge  0)
# State B: TP1  (-OPO3H-,  charge -1)
#
# No RDKit, no ITP preparation, no FF patching required.
# For TP2 (dianionic, charge -2): replace P1 with P2 throughout.
#
# Prerequisites:
#   Run:  bash input/fetch_input.sh   to download 1AOT and extract chain A.
# =============================================================================
set -euo pipefail

FF=charmm36m-mut
WATER=tip3p
INPUT=input/1AOT_A.pdb

# Target Tyr residue in pmx 1-based sequential numbering.
# Run `bash input/fetch_input.sh` to see the available tyrosines.
MUT_RESID=63       # adjust to your target Tyr
MUT_CODE=P1        # P1 = monoanionic TP1 (Y2P1); use P2 for dianionic TP2 (Y2P2)

# ── 0. Fetch input if not present ──────────────────────────────────────────
if [[ ! -f "$INPUT" ]]; then
    echo ">>> Fetching 1AOT from RCSB ..."
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
# Mutation code P1 selects the Y2P1 hybrid (TYR -> monoanionic phosphoTyr).
# P2 would select Y2P2 (TYR -> dianionic phosphoTyr).
echo ">>> pmx mutate: Tyr${MUT_RESID} -> TP1 ..."
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
echo "State A: Tyr${MUT_RESID}  (neutral,      charge  0)"
echo "State B: TP1${MUT_RESID}  (monoanionic,  charge -1)"
echo ""
echo "Note: add 1 K+ counter-ion for the state B simulation leg."
