#!/usr/bin/env bash
# Download a structure containing a serine to phosphorylate.
# Here we use 2O88 (c-Src SH2 domain + phosphopeptide).
# We extract chain A (SH2 domain) and rename residues as needed.
set -euo pipefail
PDB=2O88
URL="https://files.rcsb.org/download/${PDB}.pdb"

echo "Downloading ${PDB}.pdb ..."
curl -sSL "$URL" -o "${PDB}_raw.pdb"

# Extract chain A (SH2 domain) ATOM records only
awk '/^ATOM/ && substr($0,22,1)=="A" { print }
     /^TER/  && substr($0,22,1)=="A" { print; exit }' \
     "${PDB}_raw.pdb" > protein.pdb
echo "END" >> protein.pdb

rm -f "${PDB}_raw.pdb"
echo "Written: protein.pdb ($(grep -c '^ATOM' protein.pdb) ATOM records)"
echo ""
echo "Check which Ser residue you want to phosphorylate:"
grep "^ATOM" protein.pdb | awk '$4=="SER" {print $4, $6}' | sort -u
