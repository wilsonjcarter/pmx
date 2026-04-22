#!/usr/bin/env bash
# =============================================================================
# Example 7 [WIP] — De-novo NSAA: Ser → phosphoSer (SEP) via ITP→RTP patching
#
# System  : Protein with a target Ser
# Mutation: Ser -> SEP (O-phosphoserine, dianionic, -2)
# FF      : charmm36m-mut
# NSAA workflow: prepare_nsaa_ff -> mutate_nsaa -> pdb2gmx -> gentop
#
# State A: Ser  (-OH,    charge  0)
# State B: SEP  (-OPO3-, charge -2)
#
# Prerequisites:
#   - input/protein.pdb    (protein structure, e.g. from fetch_input.sh)
#   - input/sep.pdb        (SEP fragment with CH3 cap at CA)
#   - input/sep.itp        (GAFF2/CGenFF ITP for SEP)
#   See input/README_sep_params.md for how to prepare the ITP.
#
# For a simpler phosphorylation calculation using the built-in SP1 residue,
# see examples/06_phosphorylation_SP1/run.sh (no ITP preparation required).
# =============================================================================
set -euo pipefail

FF=charmm36m-mut
WATER=tip3p
INPUT=input/protein.pdb
SEP_PDB=input/sep.pdb
SEP_ITP=input/sep.itp

# Set the residue ID and name to mutate — adjust for your structure!
MUT_RESID=42       # example: Ser42 in chain A (pmx-renumbered, 1-based)
MUT_RESNAME=SEP    # three-letter code of the target NSAA

# ── 0. Fetch input if not present ──────────────────────────────────────────
if [[ ! -f "$INPUT" ]]; then
    echo ">>> Fetching structure from RCSB ..."
    bash input/fetch_input.sh
fi

if [[ ! -f "$SEP_PDB" || ! -f "$SEP_ITP" ]]; then
    echo "ERROR: Missing $SEP_PDB or $SEP_ITP"
    echo "See input/README_sep_params.md for instructions."
    exit 1
fi

# ── 1. Set GMXLIB ──────────────────────────────────────────────────────────
eval "$(python3 -c 'from pmx.gmx import set_gmxlib; import os; set_gmxlib(); print("export GMXLIB="+os.environ["GMXLIB"])')"

FFDIR=$(python3 -c "from pmx.utils import get_ff_path; print(get_ff_path('$FF'))")
echo "Force field directory: $FFDIR"

# ── 2. Prepare NSAA force-field files ──────────────────────────────────────
# prepare_nsaa_ff:
#   - reads the CH3-capped ITP fragment (sep.itp / sep.pdb)
#   - builds a proper RTP entry for SEP using Ser backbone + SEP sidechain
#   - writes SEP.rtp into the FF directory (used by pdb2gmx)
#   - writes SEP_bonded.itp (cross-boundary CHARMM<->GAFF bonded terms)
#   - does NOT require RDKit
echo ">>> prepare_nsaa_ff: building RTP and bonded supplement ..."
python3 - <<'PYEOF'
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'src'))
from pmx.nsaa import prepare_nsaa_ff

out_rtp, out_supplement = prepare_nsaa_ff(
    itp_file     = 'input/sep.itp',
    base_resname = 'SER',        # canonical AA the NSAA is derived from
    new_resname  = 'SEP',        # three-letter name of the NSAA
    ff           = 'charmm36m-mut',
    target_charge= -2,           # net charge of SEP (phosphate at pH 7)
    atomtypes_itp= None,         # set to atomtypes ITP if FF lacks GAFF types
    out_rtp      = 'SEP.rtp',    # will also be copied to FF dir
    out_supplement= 'SEP_bonded.itp',
)
print('RTP written to         :', out_rtp)
print('Bonded supplement      :', out_supplement)
PYEOF

# ── 3. Build hybrid structure ───────────────────────────────────────────────
# mutate_nsaa:
#   - reads protein.pdb and the FF RTP (now includes SEP)
#   - builds the hybrid residue (Ser sidechain in A, SEP sidechain in B)
#   - writes the hybrid MTP file (e.g. S2SEP.mtp) for gentop
echo ">>> pmx mutate_nsaa ..."
pmx mutate_nsaa \
    -f        "$INPUT" \
    -o        mutant.pdb \
    --resid   "$MUT_RESID" \
    --resname "$MUT_RESNAME" \
    -ff       "$FF" \
    --itp     "$SEP_ITP" \
    --nsapdb  "$SEP_PDB"

# The above writes:
#   mutant.pdb          — hybrid structure
#   S2SEP.mtp           — hybrid MTP (name = <src_one_letter>2<NSAA>)
#   charmm36m-mut.ff/   — SEP.rtp + SEP.arn (atom renaming for pdb2gmx)

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
# Pass both the MTP file (tells gentop about the S2SEP hybrid) and the
# bonded supplement (provides GAFF/cross-boundary bonded parameters for
# the SEP sidechain atoms in the B state).
echo ">>> pmx gentop ..."
pmx gentop \
    -p                topol.top \
    -o                pmxtop.top \
    --extra_mtp       S2SEP.mtp \
    --supplement_bonded SEP_bonded.itp

echo ""
echo "=== Done ==="
echo "NSAA FF patch      : $FFDIR/SEP.rtp"
echo "Bonded supplement  : SEP_bonded.itp"
echo "Hybrid MTP         : S2SEP.mtp"
echo "Hybrid structure   : mutant.pdb"
echo "Hybrid topology    : pmxtop.top"
echo ""
echo "State A: Ser${MUT_RESID}  (neutral, charge  0)"
echo "State B: SEP${MUT_RESID}  (phosphate, charge -2)"
echo ""
echo "Note: the system gains -2 charge in state B."
echo "Add 2 Na+ ions for the state B leg to keep the system neutral."
