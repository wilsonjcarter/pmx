#!/usr/bin/env bash
# Fetch a small protein structure that contains a surface-exposed serine.
# 1ERT (E. coli thioredoxin) contains Ser28 and Ser107 as convenient targets.
#
# Adjust the PDB ID and chain/residue selection for your own target.
set -euo pipefail

PDB=1ERT
curl -fsSL "https://files.rcsb.org/download/${PDB}.pdb" \
    | grep -E "^(ATOM|TER|END)" \
    > protein.pdb

echo "Downloaded ${PDB}.pdb -> protein.pdb"
echo "Surface serines (potential phosphorylation sites):"
grep "^ATOM" protein.pdb | awk '$4=="SER"{print $4, $5, $6}' | sort -u
