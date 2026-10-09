# Clean mutff — plan

**Goal:** a `mutff_clean/` FF set = **pristine OG base** (upstream deGrootLab/develop) + pmx
additions as **separate topic files only**. "What did pmx add?" = `ls` the overlay files.
Base files are never edited, so a new upstream release drops in cleanly and the overlay re-layers
onto ECC / amber unchanged.

## Decisions (locked)
- **Ship `mutres.rtp/.mtp`** — the default path for standard mutations; keep it (from OG).
- **DUM_ atom types stay in `atomtypes.atp` + `ffnonbonded.itp`** (as pmx writes them) — not split.
- **CHARMM is golden**, do it first; ECC and amber are re-layers of the same overlay afterward.
- Base = **upstream/develop** (freshest OG; local develop was stale). Build in `mutff_clean/`
  (sibling dir), leave the mid-Phase-4 `mutff/` tree untouched.

## Base (from upstream/develop, untouched)
merged.rtp, mutres.{rtp,mtp,arn}, mutres_dna.*, atomtypes.atp, ffbonded.itp, ffnonbonded.itp,
cmap.itp, ions.itp, tdb/hdb/r2b, water models, forcefield.itp.

## Overlay = pmx additions (topic files, glob-/include-loaded)
From `git diff upstream/develop..covalent-baseline -- charmm36m-mut.ff`:

| file | content | source |
|---|---|---|
| `acyl.rtp` | ALY BLY FLY HIB KCR PLY SLY | covalent-baseline (as-is) |
| `nos.rtp` + `nos_bonded.itp` + `cgenff_link.itp` | NOS/SONOS/sulfenic (CN0/KN0/KS0, CSO) | covalent-baseline (as-is) |
| `mutres_term.{rtp,mtp,arn}` | terminal hybrid library | covalent-baseline (as-is) |
| `ptm_phospho.rtp` | SEP TPO PTR YP1 YP2 | **extract from merged.rtp** |
| `ptm_methyl.rtp` | MLY MLYS M3L | **extract from merged.rtp** |
| `ptm_misc.rtp` | CSO CSS NEP ASPM (+ NME stays in base merged.rtp) | extract from merged.rtp |
| `glycation.rtp` | KCX | **port from `~/projects/pmx_glycans/charmm36m.ff`** |
| `ptm_extra.rtp` | CIR HYP MLZ SNC | port from pmx_glycans |

Residues added to merged.rtp by covalent-baseline (to relocate to topic files, EXCEPT NME which
must live in the base RTP): `ASPM CSO CSS M3L MLY MLYS NEP NME PTR SEP YP1 YP2`.

## Param-coverage check (before trusting a ported residue)
For each ported/relocated PTM residue: every `[atoms]` atom type and every `[bonds]`/impropers
term must resolve in the base `atomtypes.atp` / `ffbonded.itp` (or a shipped `*_bonded.itp` like
`nos_bonded.itp`). Method: collect the residue's atom types + bonded atom-type tuples, grep the
base bonded/atomtype files, list any misses → those are the params to add (likely CGenFF, into a
`pmx_bonded.itp`).

## Validation
- **Acceptance:** the clean FF regenerates the `tests/golden` fixtures **byte-identical** to
  covalent-baseline → nothing dropped, existing showcase legs stay valid.
- **Later:** rebuild a few peptide systems (pdb2gmx → grompp, both endpoints) per PTM class and
  confirm same topology/answer as before. Small pdb2gmx+grompp smoke tests to add.

## Order
1. `charmm36m-mut` — build overlay, param-check, golden-anchor. ← **in progress; Carter assesses**
2. re-layer onto `charmm36m-ecc-mut` (ECC-scale charged side chains in the overlay only).
3. amber: split inline GAFF-junction + GSH-cap + Zn params into `pmx_bonded.itp`; build out PTM
   content (port from CHARMM, retype to amber/GAFF).
