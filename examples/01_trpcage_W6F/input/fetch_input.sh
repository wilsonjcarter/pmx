#!/usr/bin/env bash
# Download 1L2Y and extract MODEL 1, chain A only.
# Requires: curl (or wget), grep, awk

set -euo pipefail
PDB=1L2Y
URL="https://files.rcsb.org/download/${PDB}.pdb"

echo "Downloading ${PDB}.pdb ..."
curl -sSL "$URL" -o "${PDB}_raw.pdb"

# Extract MODEL 1 and chain A ATOM records only
awk '
  /^MODEL[[:space:]]+1[[:space:]]*$/ { in_model=1 }
  /^ENDMDL/                          { if (in_model) exit }
  in_model && /^ATOM/ && substr($0,22,1)=="A" { print }
' "${PDB}_raw.pdb" > ${PDB}.pdb

# Add TER and END
echo "TER" >> ${PDB}.pdb
echo "END" >> ${PDB}.pdb

rm -f "${PDB}_raw.pdb"
echo "Written: ${PDB}.pdb ($(grep -c '^ATOM' ${PDB}.pdb) ATOM records)"
