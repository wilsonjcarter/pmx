#!/usr/bin/env bash
# Download 1ERT and extract chain A ATOM records only (no HETATM).
set -euo pipefail
PDB=1ERT
URL="https://files.rcsb.org/download/${PDB}.pdb"

echo "Downloading ${PDB}.pdb ..."
curl -sSL "$URL" -o "${PDB}_raw.pdb"

awk '/^ATOM/ && substr($0,22,1)=="A" { print }
     /^TER/  && substr($0,22,1)=="A" { print; exit }' \
     "${PDB}_raw.pdb" > ${PDB}.pdb

echo "END" >> ${PDB}.pdb
rm -f "${PDB}_raw.pdb"
echo "Written: ${PDB}.pdb ($(grep -c '^ATOM' ${PDB}.pdb) ATOM records)"
