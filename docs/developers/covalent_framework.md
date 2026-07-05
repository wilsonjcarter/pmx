# The pmx covalent-modification framework

*Design document — the unified core for arbitrary covalent modifications.*

## Thesis

The entire universe of covalent modifications is expressed by **one intermediate
representation** (a dual-state hybrid residue + inter-residue bond-morphs) and
**three composable operations** on it. Every modification is a *composition* of
these — there is no per-feature engine, only a per-feature recipe.

**The three operations** (each a small transform on the shared IR):

| Operation | Meaning | Shared machinery |
|---|---|---|
| `mutate(res, target)` | morph one residue WT→target | graph-MCS map → dummy gen → emit |
| `decouple(res, atoms)` | split reacting atoms into two branches so a forming bond's topological exclusions are harmless (dummy-in-off-state) | the same dummy machinery |
| `bond(atomA, atomB, mode)` | toggle an inter-residue bond (+induced angles/dihedrals) A↔B | `apply_bond_morph` |

**Every modification is a composition:**

| Modification | Composition |
|---|---|
| Point mutation / PTM / NSAA | `mutate` |
| Disulfide form/break | `mutate×2 + decouple×2 + bond` |
| NOS bridge | `mutate×2 + decouple×2 + bond×2` |
| Terminal deletion | `mutate×2` (one WT→all-dummy) `+ decouple` |

`decouple` is *not* a separate builder: it is `mutate` (the WT→bound morph
through the common core) followed by splitting the reacting morphing atoms into
their A-only and B-only halves. So a cross-link is literally "a mutation, plus
decouple the reacting atoms, plus a bond" — maximal reuse of one core. The only
per-feature code is the recipe (which atoms react, which bond forms).

## Layered architecture

```
                         ┌─────────────────────────────────────────────┐
 target front ends  ───► │  TargetResidue.from_rtp / from_structure /   │
 (RTP name / PDB /       │  from_itp    →  named atoms + coords + params │
  ITP fragment)          └─────────────────────────────────────────────┘
                                          │
                                          ▼
                         ┌─────────────────────────────────────────────┐
 graph mapping      ───► │  map_residues  (2-D heavy-atom MCS,          │
 (topology-guaranteed)   │  name-prior; H by mapped parent)            │
                         └─────────────────────────────────────────────┘
                                          │
                                          ▼
                         ┌─────────────────────────────────────────────┐
 hybrid build/emit  ───► │  build_hybrid_residue  (dummy gen +          │
                         │  _make_rtp / _make_mtp)   → nsaa.rtp/.mtp    │
                         └─────────────────────────────────────────────┘
                                          │
      inter-residue bond morphs           ▼
      ┌──────────────────────────┐   ┌──────────────────────────────────┐
      │ CovalentEdit /           │──►│  gen_hybrid_top  (B-state fill +  │
      │ InterResidueBondMorph /  │   │  apply_bond_morph over [inter_    │
      │ apply_bond_morph         │   │  bonds])            → pmxtop.top   │
      └──────────────────────────┘   └──────────────────────────────────┘
```

### Status
- **DONE:** target front ends (`covalent/target.py`), graph mapper
  (`covalent/mapping.py`), hybrid build/emit + single-file `nsaa.rtp/.mtp`
  (`covalent/mutate.py`), inter-residue IR + applier (`covalent/edit.py`,
  verified byte-equivalent to `alchemy._add_ss_bonded_terms`).
- **DONE:** the generic graph-derived bond builder `covalent.edit.
  induced_bond_morph` (Level 2 — declare the bond, derive angles/dihedrals from
  the union graph; replaces the hand-authored `disulfide_morph`/`nos_morph`).
- **DONE:** the `[ inter_bonds ]` MTP section — writer (`nsaa._make_mtp`),
  parser (`mutdb.read_new_mtp_entry` → `mol.inter_bonds`), generation
  (`CovalentEdit.emit`), and the additive `gen_hybrid_top` consumer
  (`alchemy._apply_inter_bonds`).  Shipped hybrids without the section are
  untouched (zero regression).
- **NEXT:** ΔG-validate the consumer against the split-atom `C2CD` reference
  (the 4 H-torsions), *then* retire `_add_ss_bonded_terms` / `_add_nos_terms` /
  `_terminal_deletion_decouplings`; multi-cross-link partner pairing.
- **PLANNED:** bond-forming exclusion handling (below).

## The exclusion constraint (why bond-forming is special)

GROMACS generates non-bonded **exclusions and 1-4 pairs topologically** from the
bond graph — *state-independently*. Forming an inter-residue bond (a disulfide
`SG–SG'`, a NOS `NZ–O`) therefore excludes the new 1-4 neighbours in **both**
alchemical states, even if the bond's force constant is 0 in state A.

Consequence: a naive "map `CYS→CYS2` and add the `SG–SG` bond" is **physically
wrong** — in state A (free thiols, no disulfide) the real `SG` loses its true
1-4 interactions with the partner residue. A single shared `SG` cannot be both
"bonded (excluded) in B" and "free (interacting) in A".

Two correct solutions:

### (a) Split-atom (current, shipped prebuilt)
Duplicate the reacting sidechain into two branches. In `C2CD`:
`CB1/HB11/HB12 + SG1/HG1` are the state-A free-thiol branch (dummy in B);
`CB2/HB21/HB22 + SG2` are the state-B disulfide branch (dummy in A). The
`SG2–SG2'` bond sits on the branch that is a **dummy in state A**, so its
topological exclusions are harmless there. Required today; **not** producible by
a 1:1 graph map (it duplicates atoms). These hybrids are shipped prebuilt (see
`disulfide.mtp`, routed via `utils.get_mtp_file` for `C2CD/C2VD/…`).

### (b) Vsite + [ pairs ] (chosen direction — lighter, generatable)
Instead of duplicating the whole sidechain, per reacting atom (`SG`, and `CB`
as needed):

- the **free-state charge** is carried by a **virtual site** co-located with
  the atom but **not** part of the SG–SG bond graph (so it is not excluded);
- the **real atom** carries the **bound-state charge** (charge morphs A→B);
- the LJ interactions removed by the bond's topological exclusions are **added
  back via explicit `[ pairs ]`** entries (state-resolved).

Advantages: only a vsite + a few pairs per reacting atom (not a duplicated
sidechain), and it can be **generated from the inter-residue bond definition**
("detect vsites/pairs from the interbond") rather than hand-built — pushing more
of disulfide/NOS into the shared core.

Open details to work out: exact vsite construction (position/type) and the
precise `[ pairs ]` set + their A/B parameters; how the vsite + pairs are
serialized alongside `[ inter_bonds ]`; validation that the resulting ΔG matches
the split-atom reference.

## Inter-residue bond serialization: `[ inter_bonds ]`

**Decision (2026-07-05, supersedes the earlier "deferred"):** the hybrid MTP
carries an `[ inter_bonds ]` section, and `gen_hybrid_top` applies it internally.
The earlier objection — that the bond references specific `Atom` *objects* lost
when the CLI reparses the topology — is dissolved by serializing the bond **by
atom name**: `gen_hybrid_top` re-resolves the names against the freshly parsed
topology and applies the morph, exactly as the resname detectors already do,
except the *"which atoms + which mode"* knowledge now lives in the file. This
removes the hardcoded resname lists and per-chemistry atom-name knowledge, so
**any** cross-link works through the CLI with no bespoke detector.

### Level 2: declare only the bond; derive the rest from the graph

The section declares **nothing but the forming bond(s) and a mode**. The induced
angles and dihedrals are generated by the **same graph rule GROMACS uses** —
every angle `x–a–b`/`a–b–y` incident on a new edge and every proper dihedral
whose path traverses it — emitted in dual (A/B) state form by the generic
builder `covalent.edit.induced_bond_morph`. This *replaces* the hand-authored
`disulfide_morph` / `nos_morph`.

```
[ inter_bonds ]
; mode   self-atom   partner-atom
  form   SG2         @SG2          ; disulfide: one bond (both Cys declare it)
```

A NOS bridge is likewise a single line (`form NZ2 @OD2`); the entire NOS
interaction set falls out of the graph from that one edge.

**Verified** (`tests/test_covalent_induced.py`):

- **NOS** — on the real `C2CN`/`K2KN` connectivity the graph reproduces
  `nos_morph` *exactly* (3 angles, 6 dihedrals): that hand-list was complete.
- **Disulfide** — the graph reproduces `disulfide_morph`'s angles exactly and its
  dihedrals *plus* the `HB-CB-SG-SG'` torsions the hand-list omitted. These are
  not zero-parameter: CHARMM36m gives `SM SM CT2 HA2` a real kφ = 0.661 kJ/mol
  (vs. 1.297 for the included `SM SM CT2 CT1`), so the shipped `C2CD` state-B
  endpoint did not exactly match a natively-built disulfide. The graph rule
  fixes this.
- **SONOS** — a two-bond bridge (`NZ–ODa`, `NZ–ODb` meeting at the Lys N): the
  graph generates the shared-nitrogen `ODa–NZ–ODb` cross-angle a per-bond
  hand-list would miss.
- **Covalent inhibitors (e.g. a BTK acrylamide → Cys thioether)** — a ligand
  atom is just another partner; the forming `SG–Cβ` bond generates every bonded
  term across the protein–ligand boundary with no per-chemistry code.

**Physics gate:** because retiring the detectors would make the shipped `C2CD`
disulfide gain those 4 H-torsions, the `gen_hybrid_top` consumer is added
**additively**; the three bespoke `alchemy.py` routines are retired only after
the resulting ΔG is validated byte-for-byte (modulo the H-torsions) against the
split-atom reference.

Force constants are resolved from the force field at apply time; the section
declares only participants + mode (`form`/`break`/`perma`/`decouple`).

## What is shared vs feature-specific

| Layer | Shared across all covalent mods? |
|---|---|
| Target front ends + graph mapper + hybrid emit | ✅ |
| Inter-residue bond formalism (`apply_bond_morph`, `[ inter_bonds ]`, consumer) | ✅ (disulfide / NOS / terminal) |
| Split-atom disulfide/NOS *residues* | ❌ physics-required; shipped prebuilt |
| Vsite + `[ pairs ]` generation (if adopted) | ✅ generatable from the bond morph |

Terminal deletion forms **no** new bond, so it has no exclusion issue and folds
into the core most directly (target → all-dummy, neighbour → cap dummies, plus a
`decouple` morph over boundary dihedrals).

## Roadmap

1. `[ inter_bonds ]` MTP section: format + parser + writer.
2. `gen_hybrid_top` consumer applying `apply_bond_morph`; retire the three
   bespoke `alchemy.py` functions (additive → validate byte-equivalence →
   remove).
3. Terminal deletion through the core (no exclusions).
4. Vsite + `[ pairs ]` generation for bond-forming (disulfide/NOS); validate ΔG
   against the split-atom prebuilt reference. Keep split-atom prebuilt as the
   fallback/reference.
5. End-to-end examples per modification class; retire superseded code.
