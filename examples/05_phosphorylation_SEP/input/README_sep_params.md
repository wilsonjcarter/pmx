# Preparing SEP (phosphoserine) parameters

You need two files:

## 1. `sep.pdb` — SEP structure with CH₃ cap at Cα

The ITP fragment must represent SEP as an isolated residue **capped** at the
backbone Cα with a methyl group (replacing the backbone C–N chain).

Atom naming convention (Cα cap):
```
  CH3(cap)–CA–CB–OG–P(=O)(OH)₂
```

Suggested preparation:
```bash
# Option A: from RCSB CCD
# Download the SEP SDF from https://www.rcsb.org/ligand/SEP
# Open in Avogadro or Maestro, add CH3 at the alpha position, minimize

# Option B: Open Babel from SMILES
#   SMILES for SEP (N-capped at CA with methyl):
#   CC(COP(=O)(O)O)N   <- too simplified; use a proper capped fragment
obabel -:"CC([NH3+])COP(=O)([O-])[O-]" --gen3d -O sep.pdb
# Then manually adjust the cap to CH3 (remove NH3+)
```

The key requirement: **the cap carbon (usually labeled C or CAP) must be
identifiable** as the CH₃ group from the ITP connectivity.

## 2. `sep.itp` — GAFF2 or CGenFF ITP for SEP

### Using GAFF2 (via Antechamber + ACPYPE)
```bash
# 1. Protonate and optimize with AM1-BCC charges
antechamber -i sep.pdb -fi pdb -o sep.mol2 -fo mol2 \
            -c bcc -nc -2 -at gaff2

# 2. Generate prmtop/inpcrd
parmchk2 -i sep.mol2 -f mol2 -o sep.frcmod -s 2

# 3. Convert to GROMACS ITP
acpype -i sep.mol2 -c bcc -n -2 -a gaff2
# -> produces sep_GMX.itp (rename to sep.itp)
```

### Using CGenFF (via ParamChem)
Upload sep.pdb to https://cgenff.umaryland.edu/ and download the `.str` file.
Then convert with `cgenff_charmm2gmx.py` or equivalent.

### Checking the ITP
Make sure:
- Atom names are unique within the ITP (no two atoms named "H1")
- The cap CH₃ atoms are present in the `[ atoms ]` section
- `[ bonds ]`, `[ angles ]`, `[ dihedrals ]` sections are present

## Automatic sanitization

If your ITP uses abstract atom types (e.g., `AT_1`, `AT_2` from some tools),
`prepare_nsaa_ff` will call `sanitize_itp` automatically to normalize names
and map types to GAFF equivalents using σ/ε nearest-neighbour matching.
