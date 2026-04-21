# Example 2 — Thioredoxin CYS→CYM protonation (1ERT)

## System

**PDB**: 1ERT — E. coli thioredoxin, reduced form (108 residues, chain A)  
**Mutation**: Cys32 → Cym  (pmx hybrid residue `C2CM`)

Thioredoxin has an unusual active-site cysteine (Cys32) with a depressed pKa
(~6.3 vs normal ~8.3) due to the surrounding electrostatic environment.
This example computes the free energy cost of deprotonating Cys32 (CYS→CYM)
inside the folded protein vs. in solution, giving the ΔpKa shift.

**CYS** = protonated thiol (-SH, neutral)  
**CYM** = deprotonated thiolate (-S⁻, charge −1)  
Hybrid residue `C2CM`: state A = CYS, state B = CYM

## Force field

`amber99sb-star-ildn-mut`

## Files

```
input/
    1ERT.pdb        # chain A of 1ERT, no HETATM
    fetch_input.sh  # download & strip
```

## Residue numbering note

1ERT has no insertion codes and residues are numbered 1–108 in the PDB file.
After pmx renumbering (default), Cys32 is still residue 32.

## Expected output

- `mutant.pdb`: hybrid structure — Cys32 sidechain unchanged in state A;
  state B is identical but carries CYM atom types (SH→S⁻, HG removed)
- `pmxtop.top`: charges morph from CYS to CYM (total charge shifts by −1)

## Notes

- The overall system becomes charged in state B; make sure to add counter-ions
  (K⁺ or Na⁺) relative to the target protonation state before production.
- Run this calculation both in explicit water (reference/alchemical leg) and
  in the folded protein to get ΔΔG = ΔpKa shift × RT ln(10).
- Also consider Cys35 at the same active site — it has a near-normal pKa (~8).
