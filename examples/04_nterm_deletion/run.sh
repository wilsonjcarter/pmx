#!/usr/bin/env bash
# =============================================================================
# Example 4 — N-terminal deletion: MEEVD peptide (1ELR chain B)
#
# System  : HSP90 C-terminal MEEVD pentapeptide
# Mutation: delete Met1 from N-terminus  ->  EdeN hybrid placed on Glu2
# FF      : charmm36m-mut
#
# State A: M(1)-E(2)-E(3)-V(4)-D(5)   (full peptide, Met1 present)
# State B:     [E(2)-E(3)-V(4)-D(5)]  (truncated peptide, Met1 = dummies)
#
# Call pmx mutate on Met1 (--resname DEL).  The code identifies Met1 as the
# N-terminus and internally renames its neighbor Glu2 to EdeN, adding dummy
# protons for the new N-terminus in state B.  Met1 stays physically present
# and is fully dummified later by pmx gentop.
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

FFDIR=$(python3 -c "from pmx.utils import get_ff_path; print(get_ff_path('$FF'))")
echo "Force field directory: $FFDIR"

# ── 2. Generate terminal-deletion parameters (once per FF) ─────────────────
# This creates mutres_term.rtp and mutres_term.mtp in the current directory.
# If write access to FFDIR is available, --outdir "$FFDIR" instead.
if [[ ! -f mutres_term.mtp ]]; then
    echo ">>> Generating terminal-deletion RTP/MTP entries ..."
    python3 -m pmx.scripts.generate_term_deletion \
        --ffdir  "$FFDIR" \
        --outdir .
fi

# Copy the generated files into the FF directory so pdb2gmx can find the RTP.
# (Alternatively, pass -r mutres_term.rtp to pdb2gmx if your version supports it.)
cp mutres_term.rtp "$FFDIR/mutres_term.rtp"
cp mutres_term.mtp "$FFDIR/mutres_term.mtp"
echo "Copied mutres_term.rtp / .mtp into $FFDIR"

# ── 3. Build hybrid structure ───────────────────────────────────────────────
# Target the residue TO BE DELETED (Met1) with --resname DEL.
# pmx mutate detects that Met1 is the N-terminus and calls
# apply_terminal_deletion(), which:
#   - finds the next residue (Glu2) and renames it to EdeN
#   - adds two dummy NH3 protons (DH2, DH3) to Glu2
#   - leaves Met1 physically present (it is fully dummified later by gentop)
echo ">>> pmx mutate: Met1 -> DEL (N-terminal deletion) ..."
pmx mutate \
    -f  "$INPUT" \
    -o  mutant.pdb \
    --resid   1 \
    --resname DEL \
    -ff "$FF"

# ── 4. Generate standard GROMACS topology ──────────────────────────────────
echo ">>> gmx pdb2gmx ..."
gmx pdb2gmx \
    -f      mutant.pdb \
    -o      processed.gro \
    -p      topol.top \
    -ff     "$FF" \
    -water  "$WATER" \
    -ignh

# ── 5. Fill hybrid B states ─────────────────────────────────────────────────
# Pass the term-deletion MTP file so gentop can resolve the EdeN residue.
echo ">>> pmx gentop ..."
pmx gentop \
    -p          topol.top \
    -o          pmxtop.top \
    --extra_mtp mutres_term.mtp

echo ""
echo "=== Done ==="
echo "Hybrid structure : mutant.pdb"
echo "Hybrid topology  : pmxtop.top"
echo ""
echo "State A: full MEEVD pentapeptide"
echo "State B: truncated EEVD tetrapeptide (Met1 = non-interacting dummies)"
