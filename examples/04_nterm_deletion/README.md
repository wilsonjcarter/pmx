# Example 4 — N-terminal residue deletion: HSP90 MEEVD peptide (1ELR)

## System

**PDB**: 1ELR — Hop TPR2A domain (chain A) + HSP90 C-terminal MEEVD peptide (chain B)  
**Mutation**: delete Met1 from the MEEVD peptide

The last 5 residues of HSP90 (Met-Glu-Glu-Val-Asp, MEEVD) bind the TPR2A
domain of Hop.  This example computes the free energy cost of truncating the
peptide from the N-terminus (removing Met1) relative to the full-length peptide.

**Terminal deletion** is a special pmx mutation that alchemically converts a
peptide into its N-terminally truncated form:

- State A: full peptide  M(1)-E(2)-E(3)-V(4)-D(5)
- State B: truncated peptide  [M(1) deleted]  E(2)-E(3)-V(4)-D(5)

The hybrid residue is placed on the **residue that becomes the new N-terminus
in state B** — here Glu2 (`EdeN`).  Met1 itself is "dummified" in state B
(all atoms present but uncharged and non-interacting).

## Force field

`charmm36m-mut` — used here to demonstrate the CHARMM workflow; the terminal
deletion MTP/RTP must be pre-generated for this force field.

## Prerequisites: generate terminal-deletion parameters

The `XdeN` hybrid residue entries are **not** built into the standard mutres.mtp;
they must be generated once per force field using:

```bash
python3 -m pmx.scripts.generate_term_deletion \
    --ffdir $(python3 -c "from pmx.utils import get_ff_path; print(get_ff_path('charmm36m-mut'))") \
    --outdir .
```

This writes `mutres_term.rtp` and `mutres_term.mtp` into the current directory
(or directly into the FF directory if you have write access).

## Files

```
input/
    1ELR_peptide.pdb   # chain B (MEEVD peptide) extracted from 1ELR
    fetch_input.sh
```

## Residue numbering

1ELR chain B runs Met401–Asp405 in PDB numbering.  After pmx renumbering
(default, starting from 1) they become 1–5.

## Expected output

- `mutant.pdb`: Glu2 (residue 2 of chain B) bears the `EdeN` hybrid type;
  Met1 is left in the structure (will be dummified at the topology level).
- `pmxtop.top`: Met1 atoms have `DUM_` prefix types in state B; Glu2 N-term
  patch charges applied in state B.
