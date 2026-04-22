#!/usr/bin/env bash
# =============================================================================
# Example 5 — Tyrosine phosphorylation: Tyr → phosphoTyr (YP1, monoanionic)
#
# System  : Lck SH2 domain (PDB 1AOT, chain A)
# Mutation: Tyr -> YP1 (monoanionic O-phosphotyrosine, charge -1)
# FF      : charmm36m-mut  (YP1 is already parameterised in this FF)
#
# State A: Tyr  (-OH,      charge  0)
# State B: YP1  (-OPO3H-,  charge -1)
#
# IMPORTANT — charge-changing mutation:
#   This transformation shifts the box charge by -1.  To avoid finite-size
#   PBC artefacts the single-box double-system (doublebox) approach is used:
#   the protein leg (Tyr→YP1, Δq=-1) and a reference peptide leg (YP1→Tyr,
#   Δq=+1) run simultaneously in the same box so the net charge never changes.
#   See README.md §Step 5 and the REF_* block below.
#
# No RDKit, no ITP preparation, no FF patching required.
# For YP2 (dianionic, charge -2): replace P1 with P2 throughout.
#
# Prerequisites:
#   Run:  bash input/fetch_input.sh   to download 1AOT, extract chain A, and
#         extract the bound Tyr-peptide as input/1AOT_peptide.pdb.
#   The reference is the BIOLOGICAL peptide ligand, not a minimal ACE-Tyr-NME.
# =============================================================================
set -euo pipefail

FF=charmm36m-mut
WATER=tip3p
INPUT=input/1AOT_A.pdb
REF_INPUT=input/1AOT_peptide.pdb  # Tyr-peptide extracted from 1AOT (biological SH2 ligand, free in water)

# Target Tyr residue in pmx 1-based sequential numbering.
# Run `bash input/fetch_input.sh` to see the available tyrosines.
MUT_RESID=63       # adjust to your target Tyr
MUT_CODE=P1        # P1 = monoanionic YP1 (Y2P1); use P2 for dianionic YP2 (Y2P2)
REF_RESID=1        # residue number of Tyr in the reference peptide

# ── 0. Fetch input if not present ──────────────────────────────────────────
if [[ ! -f "$INPUT" ]]; then
    echo ">>> Fetching 1AOT from RCSB ..."
    bash input/fetch_input.sh
fi

# ── 1. Set GMXLIB ──────────────────────────────────────────────────────────
eval "$(python3 -c 'from pmx.gmx import set_gmxlib; import os; set_gmxlib(); print("export GMXLIB="+os.environ["GMXLIB"])')"

# ════════════════════════════════════════════════════════════════════════════
# PROTEIN LEG — Y2P1 hybrid in the SH2 domain
# ════════════════════════════════════════════════════════════════════════════

# ── 2. Normalise atom names with pdb2gmx ───────────────────────────────────
echo ">>> gmx pdb2gmx (protein wildtype) ..."
gmx pdb2gmx \
    -f      "$INPUT" \
    -o      wt.gro \
    -p      wt.top \
    -ff     "$FF" \
    -water  "$WATER" \
    -ignh

# ── 3. Build hybrid structure ───────────────────────────────────────────────
echo ">>> pmx mutate: Tyr${MUT_RESID} -> YP1 ..."
printf "%s %s\n" "$MUT_RESID" "$MUT_CODE" | pmx mutate \
    -f      wt.gro \
    -o      mutant.pdb \
    -ff     "$FF"
rm -f wt.gro wt.top

# ── 4. Hybrid topology ──────────────────────────────────────────────────────
echo ">>> gmx pdb2gmx (protein mutant) ..."
gmx pdb2gmx \
    -f      mutant.pdb \
    -o      processed.gro \
    -p      topol.top \
    -ff     "$FF" \
    -water  "$WATER"

echo ">>> pmx gentop (protein) ..."
pmx gentop \
    -p  topol.top \
    -o  pmxtop.top \
    -ff "$FF"

# ════════════════════════════════════════════════════════════════════════════
# REFERENCE LEG — same Y2P1 hybrid applied to the free Tyr-peptide in water
# The reference is the BIOLOGICAL SH2 ligand (extracted from 1AOT) in free
# solution — not a minimal ACE-Tyr-NME.  Running the reverse transformation
# (YP1→Tyr, Δq=+1) keeps the combined box charge constant, and the combined
# work directly yields ΔΔG of binding: how much more tightly SH2 recruits the
# phosphopeptide vs the unphosphorylated peptide.
# ════════════════════════════════════════════════════════════════════════════

echo ">>> gmx pdb2gmx (Tyr-peptide, free in water) ..."
gmx pdb2gmx \
    -f      "$REF_INPUT" \
    -o      ref_wt.gro \
    -p      ref_wt.top \
    -ff     "$FF" \
    -water  "$WATER" \
    -ignh

echo ">>> pmx mutate: peptide Tyr${REF_RESID} -> YP1 ..."
printf "%s %s\n" "$REF_RESID" "$MUT_CODE" | pmx mutate \
    -f      ref_wt.gro \
    -o      ref_mutant.pdb \
    -ff     "$FF"
rm -f ref_wt.gro ref_wt.top

echo ">>> gmx pdb2gmx (reference mutant) ..."
gmx pdb2gmx \
    -f      ref_mutant.pdb \
    -o      ref_processed.gro \
    -p      ref_topol.top \
    -ff     "$FF" \
    -water  "$WATER"

echo ">>> pmx gentop (reference) ..."
pmx gentop \
    -p  ref_topol.top \
    -o  ref_pmxtop.top \
    -ff "$FF"

# ════════════════════════════════════════════════════════════════════════════
# DOUBLEBOX — combine protein and reference in one charge-neutral box
# ════════════════════════════════════════════════════════════════════════════

echo ">>> pmx doublebox ..."
pmx doublebox \
    -f1 processed.gro \
    -f2 ref_processed.gro \
    -o  doublebox.gro \
    -r  2.5 \
    -d  1.5

# Combine topologies: the doublebox.top includes both protein and reference
# [ molecules ] sections.  Merge pmxtop.top + ref_pmxtop.top manually or
# with a wrapper script before solvating.
# NOTE: topology merging for protein+protein doublebox is not yet automated;
#       see README.md §Step 5c for guidance.

echo ">>> gmx solvate ..."
gmx solvate \
    -cp doublebox.gro \
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
echo "Combined structure : doublebox.gro  (SH2 domain + free Tyr-peptide)"
echo "Hybrid topology    : pmxtop.top"
echo ""
echo "Protein leg (SH2 domain context)"
echo "  State A: Tyr${MUT_RESID}  (neutral, charge 0)"
echo "  State B: YP1${MUT_RESID}  (monoanionic, charge -1)"
echo "Reference leg (Tyr-peptide free in water)  [runs as B→A: YP1→Tyr]"
echo "  State A: Tyr  (neutral, charge 0)"
echo "  State B: YP1  (monoanionic, charge -1)"
echo ""
echo "Net box charge change during NEQ transitions: 0"
echo "Combined ΔΔG = ΔG(phosphorylation in SH2 context) − ΔG(phosphorylation of free peptide)"
echo "             = ΔΔG of SH2 binding: pTyr vs Tyr peptide"
