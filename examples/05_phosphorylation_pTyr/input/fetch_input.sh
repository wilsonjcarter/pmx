#!/usr/bin/env bash
# Fetch 1AOT (Lck SH2 domain + bound phosphopeptide) and write the two inputs
# example 5 needs.
#
# 1AOT chains: F = SH2 domain (106 residues), P = ligand peptide EPQEEIPIYL.
# There is NO chain A.  Residue/chain fields are read by PDB COLUMN, not by
# whitespace field, because a PDB record is fixed-format: altLoc, four-character
# atom names and long residue names all shift the split-field positions.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # write beside THIS script
PDB=1AOT
RAW="$HERE/${PDB}.pdb"
COMPLEX="$HERE/${PDB}_complex.pdb"
PEPTIDE="$HERE/${PDB}_peptide.pdb"

curl -fsSL "https://files.rcsb.org/download/${PDB}.pdb" -o "$RAW"

# complex = both chains, the binding-competent system
awk '/^ATOM/ && (substr($0,22,1)=="F" || substr($0,22,1)=="P")' "$RAW" > "$COMPLEX"
printf 'TER\nEND\n' >> "$COMPLEX"

# reference = the ligand peptide alone, free in water
awk '/^ATOM/ && substr($0,22,1)=="P"' "$RAW" > "$PEPTIDE"
printf 'TER\nEND\n' >> "$PEPTIDE"

echo "Wrote $(basename "$COMPLEX") (chains F+P) and $(basename "$PEPTIDE") (chain P)"
echo ""
echo "Tyrosines in the ligand peptide, in pmx 1-based sequential numbering:"
awk '/^ATOM/ && substr($0,13,4)==" CA " {n++; if (substr($0,18,3)=="TYR") printf "  peptide residue %d (TYR)\n", n}' "$PEPTIDE"
echo ""
echo "Same tyrosine inside the complex (chain F precedes chain P):"
awk '/^ATOM/ && substr($0,13,4)==" CA " {n++; if (substr($0,18,3)=="TYR" && substr($0,22,1)=="P") printf "  complex residue %d (TYR, chain P)\n", n}' "$COMPLEX"
