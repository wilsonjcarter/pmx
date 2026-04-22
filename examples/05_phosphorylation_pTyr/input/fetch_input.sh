#!/usr/bin/env bash
# Fetch the Lck SH2 domain structure (1AOT) from the RCSB and extract chain A.
# Chain A is the SH2 domain; other chains (if present) are ligand/peptide.
set -euo pipefail

PDB=1AOT
OUTFILE=1AOT_A.pdb

curl -fsSL "https://files.rcsb.org/download/${PDB}.pdb" \
    | awk '/^ATOM/ && $5=="A" {print} /^TER/ {print; exit}' \
    > "$OUTFILE"

echo "Downloaded ${PDB}.pdb -> ${OUTFILE} (chain A only)"
echo ""
echo "Tyrosine residues available as phosphorylation targets:"
grep "^ATOM" "$OUTFILE" | awk '$4=="TYR" && $3=="CA" {print "  TYR", $6, "(chain", $5")"}' | sort -u
echo ""
echo "Use the residue number (column 6) in pmx 1-based sequential numbering."
echo "Run: awk '/^ATOM/{print NR, \$4, \$6}' $OUTFILE | grep TYR | head -20"
echo "to see 1-based atom-record positions as a cross-check."
