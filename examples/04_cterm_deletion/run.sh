#!/usr/bin/env bash
# =============================================================================
# Example 4 — C-terminal deletion: MEEVD peptide (1ELR chain B)
#
# System  : HSP90 C-terminal MEEVD pentapeptide
# Mutation: delete Asp5 from C-terminus  ->  VdeC hybrid placed on Val4
# FF      : charmm36m-mut
#
# State A: M(1)-E(2)-E(3)-V(4)-D(5)   (full peptide, Asp5 present)
# State B: M(1)-E(2)-E(3)-V(4)         (truncated peptide, Asp5 = dummies)
#
# Call pmx mutate on Asp5 (--resname DEL).  The code identifies Asp5 as the
# C-terminus and internally renames its neighbour Val4 to VdeC, adding a
# dummy oxygen (DOT2) that represents OT2 of the new C-terminal carboxylate
# in state B.  Asp5 stays physically present and is fully dummified later
# by pmx gentop.
# =============================================================================
set -euo pipefail

FF=charmm36m-mut
WATER=tip3p
INPUT=input/1ELR_peptide.pdb

# ── 0. Fetch input if not present ──────────────────────────────────────────
if [[ ! -f "$INPUT" ]]; then
    echo ">>> Fetching 1ELR peptide chain from RCSB ..."
    bash input/fetch_input.sh
fi

# ── 1. Set GMXLIB ──────────────────────────────────────────────────────────
eval "$(python3 -c 'from pmx.gmx import set_gmxlib; import os; set_gmxlib(); print("export GMXLIB="+os.environ["GMXLIB"])')"

FFDIR=$(python3 -c "import os; from pmx.gmx import set_gmxlib; set_gmxlib(); print(os.environ['GMXLIB']+'/charmm36m-mut.ff')")
echo "Force field directory: $FFDIR"

# ── 2. Generate terminal-deletion parameters (once per FF) ─────────────────
# This creates mutres_term.rtp and mutres_term.mtp in the current directory.
if [[ ! -f mutres_term.mtp ]]; then
    echo ">>> Generating terminal-deletion RTP/MTP entries ..."
    python3 -m pmx.scripts.generate_term_deletion \
        --ffdir  "$FFDIR" \
        --outdir .
fi

# Copy the RTP into the FF directory so pdb2gmx can find VdeC.
cp mutres_term.rtp "$FFDIR/mutres_term.rtp"
cp mutres_term.mtp "$FFDIR/mutres_term.mtp"
echo "Copied mutres_term.rtp / .mtp into $FFDIR"

# ── 3. Normalise atom names with pdb2gmx ───────────────────────────────────
echo ">>> gmx pdb2gmx (wildtype — atom-name normalisation) ..."
gmx pdb2gmx \
    -f      "$INPUT" \
    -o      wt.gro \
    -p      wt.top \
    -ff     "$FF" \
    -water  "$WATER" \
    -ignh

# ── 4. Build hybrid structure ───────────────────────────────────────────────
# Target the residue TO BE DELETED (Asp5) with "5 DEL".
# pmx mutate detects that Asp5 is the C-terminus and calls
# apply_terminal_deletion(), which:
#   - finds the previous residue (Val4) and renames it to VdeC
#   - adds a dummy oxygen (DOT2) to Val4 for the future OT2 C-terminal carboxylate
#   - leaves Asp5 physically present (it is fully dummified later by gentop)
echo ">>> pmx mutate: Asp5 -> DEL (C-terminal deletion) ..."
printf "5 DEL\n" > mut.txt
pmx mutate \
    -f      wt.gro \
    -o      mutant.pdb \
    --script mut.txt \
    -ff     "$FF"
rm -f mut.txt wt.gro wt.top

# ── 5. Generate GROMACS topology for the hybrid structure ──────────────────
# mutres_term.rtp has been copied into FFDIR so pdb2gmx can find VdeC.
# Do NOT pass -ignh: the dummy atom DOT2 placed by pmx mutate must be
# preserved; pdb2gmx would not know how to re-add it from HDB.
echo ">>> gmx pdb2gmx (mutant) ..."
gmx pdb2gmx \
    -f      mutant.pdb \
    -o      processed.gro \
    -p      topol.top \
    -ff     "$FF" \
    -water  "$WATER"

# ── 6. Fill hybrid B states ─────────────────────────────────────────────────
# Pass the term-deletion MTP file so gentop can resolve the VdeC residue.
echo ">>> pmx gentop ..."
pmx gentop \
    -p          topol.top \
    -o          pmxtop.top \
    -ff         "$FF" \
    --extra_mtp mutres_term.mtp

echo ""
echo "=== Done ==="
echo "Hybrid structure : mutant.pdb"
echo "Hybrid topology  : pmxtop.top"
echo ""
echo "State A: full MEEVD pentapeptide"
echo "State B: truncated MEEV tetrapeptide (Asp5 = non-interacting dummies)"
