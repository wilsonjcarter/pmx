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
- **CHOSEN:** bond-forming exclusion handling via generated split-atom
  dual-topology (`decouple`); see below.

## The exclusion constraint (why bond-forming is special)

GROMACS generates non-bonded **exclusions and 1-4 pairs topologically** from the
bond graph — *state-independently*. Forming an inter-residue bond (a disulfide
`SG–SG'`, a NOS `NZ–O`) therefore excludes the new 1-4 neighbours in **both**
alchemical states, even if the bond's force constant is 0 in state A.

Consequence: a naive "map `CYS→CYS2` and add the `SG–SG` bond" is **physically
wrong** — in state A (free thiols, no disulfide) the real `SG` loses its true
1-4 interactions with the partner residue. A single shared `SG` cannot be both
"bonded (excluded) in B" and "free (interacting) in A".

### The solution: generated split-atom dual-topology (`decouple`)

Duplicate the reacting sidechain into two decoupled branches. In `C2CD`:
`CB1/HB11/HB12 + SG1/HG1` are the state-A free-thiol branch (dummy in B);
`CB2/HB21/HB22 + SG2` are the state-B disulfide branch (dummy in A). The
`SG2–SG2'` bond sits on the branch that is a **dummy in state A**, so its
topological exclusions are harmless there — the LJ/charge the real atom loses to
those exclusions is "clawed back" by the *other* branch, which is real in that
state.

This is **not** a hand-built special case: `decouple(res, reacting_atoms)` walks
the sidechain from the bond atom, dummifies, and generates exactly this
two-branch split from the same dummy machinery `mutate` uses. So the shipped
`C2CD`/`C2VD/…` split residues are reproducible — and *any* new cross-link is
too — with no per-chemistry code. It is the general mechanism, so it is the one
we build on.

**Why not virtual sites + `[ pairs ]`.** An alternative was considered: keep a
single reacting atom, carry the free-state charge on a co-located zero-LJ virtual
site (outside the bond graph, so not excluded), morph the charge onto the real
atom in the bound state, and add back the excluded LJ via explicit state-resolved
`[ pairs ]`. It is lighter on particle count, but it is **rejected**: it needs
new, fiddly machinery (vsite construction records, a hand-tuned `[ pairs ]` set
with A/B parameters, local exclusion bookkeeping) that is harder to validate and
does not generalize. The dummy-clawback split reuses the core dummy mechanism and
covers every bond-forming case uniformly, so it is strictly more general.

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
| Split-atom disulfide/NOS *residues* (`decouple`) | ✅ generated from the reacting-atom set; the shipped prebuilts are reproduced, not special |

Terminal deletion forms **no** new bond, so it has no exclusion issue and folds
into the core most directly (target → all-dummy, neighbour → cap dummies, plus a
`decouple` morph over boundary dihedrals).

**Implemented (`covalent.delete_terminal`).** A terminal deletion is now two
residue morphs composed from the engine, with *no* bespoke machinery:

* the **neighbour** morph (internal → terminus) is an ordinary
  `build_hybrid_residue` graph map of `internal → terminal`, where the terminal
  state is `TargetResidue.from_rtp(resname, term=…)` (the same
  `nsaa._rtp_terminal` patch the NSAA path uses).  It grows the terminus cap
  (`OT2` / extra amide H's) as an appearing dummy and retypes the backbone,
  emitted as the `XdeC`/`XdeN` hybrid into `nsaa.{rtp,mtp}` like any other
  mutation.  `build_hybrid_residue` gained `term_b`/`hybrid_name` so state A and
  state B can carry different terminal patches under a chosen name.
* the **target** morph (WT → all-dummy) is built through the same engine so
  `_make_bstate_dummies` + `_write_atp_fnb` register the backbone `DUM_` types
  (e.g. `DUM_NH1`) an ordinary sidechain mutation never introduces.
  gen_hybrid_top applies the structural dummification and the boundary
  `decouple` (via `apply_bond_morph`).

This retired `scripts/generate_term_deletion.py` (the hand-authored per-force-
field `XdeC`/`XdeN` generator) and `alchemy.apply_terminal_deletion` + its
dummy-coordinate builders.  Charge bookkeeping was corrected so a wholly-
dummified residue no longer aborts `Topology.write`: `get_hybrid_qA/qB` iterate
*perturbed* residues (aligning `target_qB` with `write_atoms`), and the integer-
charge correction skips a residue with no real atom to hold the residual.

*Known limitation:* the split (`gentop --split`) legs still carry the neighbour
backbone CMAP across the deleted `+N`, which the all-state-B `qon` leg rejects
(`Unknown cmap torsion … DUM`).  The shipped example uses the dual-state
topology (no `--split`), which is CMAP-clean; decoupling the boundary CMAP for
the split legs is the remaining terminal-deletion item.

## Roadmap

1. `[ inter_bonds ]` MTP section: format + parser + writer.
2. `gen_hybrid_top` consumer applying `apply_bond_morph`; retire the three
   bespoke `alchemy.py` functions (additive → validate byte-equivalence →
   remove).
3. ~~Terminal deletion through the core (no exclusions).~~ ✅ done
   (`covalent.delete_terminal`); remaining: boundary-CMAP decouple for the
   `--split` legs.
4. `decouple`-generated split residues for bond-forming (disulfide/NOS);
   validate the generated hybrids reproduce the shipped split-atom prebuilts,
   then generate rather than ship them.
5. End-to-end examples per modification class; retire superseded code.
