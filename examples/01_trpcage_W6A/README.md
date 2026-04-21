# Example 1 — Trp Cage W6A (standard substitution)

## System

**PDB**: 1L2Y — NMR ensemble of the Trp Cage miniprotein (20 residues, single chain A)  
**Sequence**: N-L-Y-I-Q-**W**-L-K-D-G-G-P-S-S-G-R-P-P-P-S  
**Mutation**: Trp6 → Ala  (pmx hybrid residue `W2A`)

This is the textbook pmx demonstration case.  The Trp residue at position 6
stabilises the hydrophobic core of the mini-protein; the ΔΔG of W6A is large
and experimentally well characterised (~2 kcal/mol destabilising).

## Force field

`amber99sb-star-ildn-mut`  — the standard pmx protein force field.

## Files

```
input/
    1L2Y.pdb          # first NMR model (MODEL 1), chain A only
    fetch_input.sh    # download & strip to MODEL 1 / chain A
```

## Workflow summary

```
pmx mutate -f input/1L2Y.pdb -o mutant.pdb --resid 6 --resname A -ff amber99sb-star-ildn-mut
gmx pdb2gmx ...
pmx gentop  ...
```

See `run.sh` for the full, executable workflow.

## Expected output

- `mutant.pdb`: hybrid structure with W6 replaced by W2A hybrid residue (has both
  Trp side-chain atoms in state A and dummy Ala atoms in state B)
- `pmxtop.top`: hybrid topology with B-state charges and atom types filled in

## Notes

- pmx renumbers residues from 1 by default, so `--resid 6` refers to the 6th
  residue in the chain regardless of the original PDB numbering.
- 1L2Y is an NMR structure; only MODEL 1 is used here.
- The `--resname A` argument passes the **one-letter code** of the target residue.
