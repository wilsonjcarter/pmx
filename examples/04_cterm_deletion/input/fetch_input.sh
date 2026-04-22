#!/usr/bin/env bash
# Download 1ELR and extract the MEEVD peptide (chain B, residues 401-405).
# The TPR domain (chain A) is retained in 1ELR_complex.pdb for reference.
set -euo pipefail
PDB=1ELR
URL="https://files.rcsb.org/download/${PDB}.pdb"

echo "Downloading ${PDB}.pdb ..."
curl -sSL "$URL" -o "${PDB}_raw.pdb"

# Peptide only (chain B) — the MEEVD 5-mer
awk '/^ATOM/ && substr($0,22,1)=="B" { print }
     /^TER/  && substr($0,22,1)=="B" { print; exit }' \
     "${PDB}_raw.pdb" > 1ELR_peptide.pdb
echo "END" >> 1ELR_peptide.pdb

# Full complex for reference
awk '/^ATOM/ { print }
     /^TER/  { print }' \
     "${PDB}_raw.pdb" > 1ELR_complex.pdb
echo "END" >> 1ELR_complex.pdb

rm -f "${PDB}_raw.pdb"
echo "Written:"
echo "  1ELR_peptide.pdb  ($(grep -c '^ATOM' 1ELR_peptide.pdb) ATOM records)"
echo "  1ELR_complex.pdb  ($(grep -c '^ATOM' 1ELR_complex.pdb) ATOM records)"
