# Example 3 — Thioredoxin disulfide formation C2CD (1ERT)

## System

**PDB**: 1ERT — E. coli thioredoxin, **reduced** form (108 residues, chain A)  
**Mutation**: Cys32 + Cys35 → both `C2CD`

Thioredoxin cycles between reduced (dithiol) and oxidised (disulfide) forms.
The active-site disulfide Cys32–Cys35 is the key chemical event.  This example
computes the free energy of disulfide **formation** in the protein relative to
a reference peptide in solution.

**C2CD**: hybrid residue for a CYS participating in forming a disulfide bond:
- State A: reduced CYS (free thiol, -SH, SG type `S`, charge −0.23)
- State B: disulfide-bonded CYS (−S–, SG type `SM`, charge −0.08; HG1 → dummy)

Both partner Cys residues must be mutated to `C2CD` in separate `pmx mutate`
calls.  `pmx gentop` automatically adds the cross-residue SG–SG bond, angles,
and dihedrals required for the hybrid topology.

## Force field

`charmm36m-mut` — the `C2CD` / `D2DC` hybrid entries and the `CYS2` residue
(the disulfide-bonded cysteine used as the state B reference) are all now
present in the bundled `mutres.rtp`, `mutres.mtp`, and `merged.rtp`.

## Files

```
input/
    1ERT.pdb        # chain A of 1ERT (reduced form — no disulfide in structure)
    fetch_input.sh
```

## Residue numbering

1ERT residues 1–108 in PDB numbering; pmx renumbers from 1 by default →
Cys32 and Cys35 retain their original residue numbers.

## Expected output

- `mutant.pdb`: both Cys32 and Cys35 are named `C2CD` hybrid residues
- `pmxtop.top`: cross-residue SG–SG bond and flanking angle/dihedral
  terms added by gentop; HG1 dummy in state B; SG type S→SM

## Notes

- Use `-ss no` with `gmx pdb2gmx` to prevent it from trying to form a
  standard (non-alchemical) disulfide bond automatically.
- The reverse transformation (disulfide breaking) uses `D2DC` for both
  partner residues.
- The SG–SG distance in 1ERT (reduced form) is ~3.9 Å; in state B it
  should be ~2.05 Å. Allow generous equilibration.
