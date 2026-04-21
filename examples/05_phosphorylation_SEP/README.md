# Example 5 — Phosphoserine (SEP) via NSAA workflow

## System

**PDB**: 2O88 — c-Src kinase SH2 domain + phosphopeptide (pY+3 position contains Ser)  
or any structure with a **serine** you want to alchemically phosphorylate.

**Mutation**: Ser → phosphoSer (SEP, O-phosphoserine)

This example demonstrates the full `pmx mutate_nsaa` + `pmx gentop` pipeline
for non-standard amino acids (NSAAs) / post-translational modifications (PTMs).
No RDKit is required; the ITP-to-RTP framework (`prepare_nsaa_ff`) handles
parameterization directly.

## Chemical context

Phosphoserine (SEP) carries a net −2 charge at physiological pH:
- State A: Ser (−OH, neutral)
- State B: SEP (−OPO₃²⁻, charge −2)

The free energy difference gives the phosphorylation affinity and its
effect on protein stability / binding.

## Force field

`charmm36m-mut` — CHARMM-compatible parameters for both the protein backbone
and the GAFF/CGenFF sidechain of SEP.

## Required input files

In addition to the protein PDB, you need:

1. **`sep.pdb`** — SEP structure with a CH₃ cap at the Cα position  
   (the methyl group replaces the backbone; Cα is the "cap carbon").  
   Source: RCSB CCD entry SEP, or build with Avogadro/Maestro/Open Babel.

2. **`sep.itp`** — GROMACS ITP fragment for SEP, parameterised with GAFF2  
   (e.g., via Antechamber/ACPYPE) or CGenFF (via ParamChem).  
   The cap CH₃ must match the cap in `sep.pdb`.

See `input/README_sep_params.md` for details on preparing these files.

## Workflow summary

```
Step 1: prepare_nsaa_ff   →  SEP.rtp + SEP_bonded.itp  (FF patch files)
Step 2: pmx mutate_nsaa   →  mutant.pdb                (hybrid structure)
Step 3: gmx pdb2gmx       →  topol.top                 (standard topology)
Step 4: pmx gentop        →  pmxtop.top                (hybrid topology)
```

## Files

```
input/
    protein.pdb            # structure with Ser to phosphorylate
    sep.pdb                # SEP fragment with CH3 cap at CA
    sep.itp                # GAFF2 ITP for SEP (from ACPYPE / antechamber)
    fetch_input.sh         # download protein structure
    README_sep_params.md   # how to prepare sep.pdb and sep.itp
```

## Expected output

- `SEP.rtp` (in `charmm36m-mut.ff/`): new residue topology for pdb2gmx
- `SEP_bonded.itp`: bonded supplement (GAFF↔CHARMM cross-boundary terms +
  ITP-internal GAFF terms)
- `A2SEP.mtp` (or similar hybrid name): hybrid MTP for gentop
- `mutant.pdb`: hybrid structure — Ser sidechain in state A,
  phosphate sidechain in state B
- `pmxtop.top`: full hybrid topology
