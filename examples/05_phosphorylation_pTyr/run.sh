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
INPUT=input/1AOT_complex.pdb      # SH2 domain (chain F) + ligand peptide (chain P)
REF_INPUT=input/1AOT_peptide.pdb  # the same ligand peptide, free in water

# The SAME tyrosine is phosphorylated in both legs -- that is what closes the
# cycle and makes the result a ddG of binding.  It is the LIGAND's Tyr, not one
# of the SH2 domain's own: chain P is EPQEEIPIYL, so Tyr is peptide residue 9,
# which in the complex (chain F's 106 residues come first) is residue 115.
# Both numbers are in pmx 1-based sequential numbering; run
# `bash input/fetch_input.sh` to print them for yourself.
MUT_RESID=115      # ligand Tyr within the complex
MUT_CODE=P1        # P1 = monoanionic YP1 (Y2P1); use P2 for dianionic YP2 (Y2P2)
REF_RESID=9        # the same Tyr in the free peptide

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
printf "%s %s\n" "$MUT_RESID" "$MUT_CODE" > mut.txt
pmx mutate \
    -f      wt.gro \
    -o      mutant.pdb \
    --script mut.txt \
    -ff     "$FF"
rm -f mut.txt
rm -f wt.gro wt.top

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

echo ">>> pmx mutate: peptide YP1 -> Tyr${REF_RESID}, REVERSED ..."
printf "%s %s\n" "$REF_RESID" "$MUT_CODE" > ref_mut.txt
pmx mutate \
    -f      ref_wt.gro \
    -o      ref_mutant.pdb \
    --script ref_mut.txt \
    --reverse \
    -ff     "$FF"
rm -f ref_mut.txt
rm -f ref_wt.gro ref_wt.top

# ════════════════════════════════════════════════════════════════════════════
# DOUBLEBOX — complex and free peptide in one charge-neutral box
#
# ORDER MATTERS: doublebox the two MUTANT STRUCTURES, then run pdb2gmx ONCE on
# the combined file.  That is what produces a single topology covering both
# legs; building a topology per leg and trying to merge them afterwards is the
# problem this ordering avoids entirely.
# ════════════════════════════════════════════════════════════════════════════

echo ">>> pmx doublebox ..."
pmx doublebox \
    -f1 mutant.pdb \
    -f2 ref_mutant.pdb \
    -o  doublebox.pdb \
    -r  2.5 \
    -d  1.5

# ── Generate ONE topology covering both legs ────────────────────────────────
# Do NOT pass -ignh: the phosphate atoms placed by pmx mutate must be kept.
# Three chains (SH2 domain, bound peptide, free peptide) -> three terminus
# selection pairs.
echo ">>> gmx pdb2gmx (doubleboxed system) ..."
printf '0\n0\n0\n0\n0\n0\n' | gmx pdb2gmx \
    -f      doublebox.pdb \
    -o      doublebox.gro \
    -p      topol.top \
    -ff     "$FF" \
    -water  "$WATER" \
    -ter

echo ">>> pmx gentop ..."
pmx gentop \
    -p  topol.top \
    -o  pmxtop.top \
    -ff "$FF"

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
