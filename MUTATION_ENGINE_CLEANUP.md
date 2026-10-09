# Mutation / hybrid engine — cleanup & unification plan

_Drafted 2026-07 during the covalent-showcase build, after hitting the same bug in
two places. Read this top-to-bottom; the plan is self-contained._

## RESUME HERE (status as of 2026-07)

**Goal in one line:** collapse the two parallel mutation orchestrations
(`covalent.mutate` + `nsaa.mutate_nsaa`) onto one engine + one toolkit, structured as
the three stages below, **without changing behaviour** — proven by a byte-identical
golden guard.

**What's already done (in the working tree, uncommitted):**
- Behaviour fixes: heavy-atom-only templates now complete from the RTP in *both*
  paths (`covalent.target.complete_target_from_structure`); robustness/clear errors
  in `mutdb.py`, `alchemy.py`. See `git diff` on `src/pmx/{covalent,nsaa,mutdb,alchemy}`.
- **Phase 0 (golden guard): BUILT.** `tests/golden/` — 13 cases (5 point mutations,
  5 protonations, terminal amidation, 2 NSAA), regenerate byte-identical. **130 tests
  pass.** Run the tripwire after every edit (exit 0 = clean):

  ```sh
  # NOTE: GMXLIB must be ABSOLUTE (the harness chdir's into a temp dir, so a
  # relative FF path stops resolving), and use the devpmx env's python + PATH
  # (has rdkit + the `pmx`/`gmx` CLIs) with PYTHONPATH forced to THIS checkout's
  # src (else `import pmx` silently picks up the off-limits installed GitHub/pmx).
  export PATH=/Users/carter/miniforge3/envs/devpmx/bin:/opt/homebrew/bin:$PATH
  export PYTHONPATH=/Users/carter/Documents/pmx/src
  export GMXLIB=/Users/carter/Documents/pmx/src/pmx/data/mutff
  PY=/Users/carter/miniforge3/envs/devpmx/bin/python
  $PY tests/golden/generate_golden.py --check   # 0 drift
  $PY -m pytest -q                               # 130 passed
  ```
  (`gmx` = /opt/homebrew/bin/gmx; `pmx` CLI = devpmx/bin/pmx. The canonical
  `charmm36m-mut.ff` working tree carries ~100 UNTRACKED leftover hybrid files
  + `nsaa.*`; the harness now strips these from its isolated temp FF copy, so a
  dirty tree no longer breaks the generative fixtures.)
- **Phase 1, step 1: DONE** — the two completion helpers unified into one
  (`covalent/target.py`). Golden 0-drift.
- **Matrix now COMPLETE (13/13 full).** The 2 NSAA cases emit `hybrid.top`: the
  failures were *fixture-environment* pollution + a documented-usage gap, not engine
  bugs. Fixes (all test-infra / a latent bug — no engine behaviour change):
  (1) `generate_golden.py` strips runtime hybrid artifacts (`nsaa.*`, `?2??.{rtp,mtp,arn}`)
  from the isolated temp FF copy — the canonical `charmm36m-mut.ff` working tree carries
  ~100 untracked leftover hybrids that collided with the freshly-written `nsaa.rtp`
  ("`K2AL` defined in two places"); (2) it now feeds any working-dir `<hybrid>.mtp` to
  `gentop --extra_mtp` (nsaa.mutate_nsaa writes its MTP to cwd, covalent accumulates into
  the FF's `nsaa.mtp`); (3) fixed `utils.py:249` `raise(ValueError, ...)` → `raise
  ValueError(...)` (Py2 tuple-raise). Guard **hardened**: a committed fixture with no
  regenerated counterpart now counts as DRIFT (was silently skipped — would have hidden a
  broken emit during the builder refactor). Both new `hybrid.top` are byte-stable across
  runs. **130 tests pass, 0 drift.** Now the target builder is guarded by full-topology
  fixtures on both the mutate_nsaa (structure) and covalent (from_rtp) paths.

**Immediate next task:** Phase 1, step 2 — **unify the target *builder***
(`covalent._build_target_residue` vs the inline load block in `nsaa.mutate_nsaa`,
`nsaa.py:~1435`). They differ in the standard-AA branch (nsaa `build_chain` vs covalent
`from_rtp`) — go one branch at a time, golden-check each. This is Stage 1 below, and it
must absorb the `.itp` CH3-Cα-cap → residue workflow and the terminal edits.

**Also open (surfaced by the golden harness — good warm-up):** the 2 NSAA golden cases
fail downstream — `mutate_nsaa` does `raise <non-exception>`
(`TypeError: exceptions must derive from BaseException`) and the from-rtp generative
path fails pdb2gmx. Fix these and their `hybrid.top` completes the matrix.

**Guardrails:** work only in `Documents/pmx` (GitHub/pmx is off-limits); use the
**canonical** `src/pmx/data/mutff` FF for golden (NOT the showcase bundle, which is
mutated); every commit must keep 130 tests + 0 golden drift.

## Why
The mutation → hybrid → emit pipeline is spread across several large modules with
**two parallel orchestrations** of one shared toolkit. Symptom that motivated this:
a heavy-atom-only crystal template broke the build, and the fix had to be written
**twice** (once in `nsaa.py`, once in `covalent/mutate.py`) because both orchestrate
target-residue construction independently. That "fix-it-twice" hazard is the thing
worth removing.

This is **not** a rewrite. `pmx/covalent/__init__.py` already states the intended
direction: covalent is a front-end facade "grown additively" that "will progressively
host [nsaa's] pieces without changing behaviour." This doc just sequences that.

## Architecture — the whole engine is three stages

Everything the engine does is one of three stages. The current code blurs Stage 1
and Stage 2 (which atoms exist vs where they sit); separating them cleanly is the
point of the cleanup.

**Stage 1 — Target residue (geometry + parameters).** Produce a single, *complete,
parameterized* residue. Two sources, one contract:
- from `.rtp`/`.itp` — build the graph, embed 3-D (RDKit);
- from a PDB residue (already the modified residue) — read heavy atoms.
Both are **heavy-atom problems**. Hydrogens are added from connectivity.
- **Rule: add hydrogens with RDKit, from the RTP/ITP graph — never rely on
  `pdb2gmx`.** These residues have no `.hdb`, so gromacs cannot place their H's
  (observed: "add HV4 to `K2AL.hdb`"; `-ignh` → "missing atoms"). Stage 1 therefore
  **emits a complete H-included residue**, and nothing downstream may depend on
  gromacs to add them. (This is exactly what the 2026-07 heavy-atom fix does — it
  should be the rule, not a patch.)

**Stage 2 — Hybrid topology (MCS map).** Merge WT and target by maximum common
substructure of the **bond graphs** — coordinate-free *in principle*. This is also
where reactive atoms are identified and the clawback / side-chain split enforces
state-independent 1–4 exclusions.
- **Caveat (the one thing that makes full-generative non-trivial):** the graph is
  coordinate-free, but the *best* mapping is not. Symmetric groups (Phe ring,
  Leu/Val methyls, Arg guanidinium, carboxylate O's) have graph-equivalent choices
  where the correct one is min-RMSD. So there is a **one-way coordinate feed**
  Stage 1 → Stage 2: pass geometry to break MCS symmetry ties. The hardcoded
  libraries encode this choice by hand; a generative path must compute it.

**Stage 3 — Bond introduction (at gentop).** Induced angles/dihedrals from the
union bond graph for a formed/broken inter-residue bond.

### Generative vs the hardcoded libraries — decide by measuring, not forking
Architecturally, full-generative is the thesis (a point mutation and an NSAA install
are the same operation). But the choice needn't be made upfront: **the Phase-0 golden
guard *is* the experiment.** Generate the 20×20 + protonation hybrids and diff against
`mutres.rtp`/`mutres.mtp`:
- byte-identical → the library is redundant; keep it only as a regression fixture;
- diffs → they name exactly what the curated library encodes that generation doesn't
  (almost certainly the symmetry mapping above) → close that gap, or keep the library
  for those cases and generate the rest.

Either way the resolution is **one generative code path, with the library demoted
from a parallel *engine* to a *cache / oracle*.** Every line of validated
infrastructure survives as a fixture; there is only one engine to maintain; and the
golden diff makes it impossible for a generative change to silently move a published
ΔG.

## Current map (grounded, with file refs)

**Entry points (3 functions, 2 CLI commands):**
- `pmx mutate` → `scripts/mutate.py` → `alchemy.mutate` (`alchemy.py:48`) — router:
  sends PTM/NSAA/terminal/generative cases to `covalent.mutate`, else classic
  `apply_aa_mutation`.
- `pmx mutate_nsaa` → `scripts/mutate_nsaa.py` → `nsaa.mutate_nsaa` (`nsaa.py:1347`)
  — the older direct entry.

**Two orchestrations, one toolkit.** `covalent/mutate.py:41` imports ~24 private
helpers from `nsaa.py` (`_merge_molecules`, `_make_bstate_dummies`, `_make_rtp`,
`_make_mtp`, `_assign_rtp_entries`, …). So the low-level toolkit lives in `nsaa.py`;
`covalent.mutate` and `nsaa.mutate_nsaa` are two orchestrations of it. covalent's own
docstring: "the emit-path counterpart of `pmx.nsaa.mutate_nsaa`."

**Verified duplication (the real target):**
1. **Target-residue building.** `covalent._build_target_residue` (`covalent/mutate.py:297`)
   vs the target-load block inside `nsaa.mutate_nsaa` (`nsaa.py:~1425–1470`). Same
   `itp > structure > rtp` precedence. *This is where the heavy-atom bug lived twice.*
2. **Completion helpers (added this session).** `nsaa._complete_target_from_rtp`
   and `covalent._complete_target_from_structure` — near-identical Kabsch-overlay.
3. **Two orchestrators** doing the same map → hybridize → emit dance.

**To assess during Phase 0 (nuanced, do NOT assume duplicated):**
- Emission: `covalent._accumulate_hybrid` (`covalent/mutate.py:360`, covalent-only)
  vs nsaa's `_make_rtp`/`_make_mtp` + its own writing. Both write hybrid entries into
  the FF `nsaa.{rtp,mtp}`, but via different code — confirm they produce identical
  output before merging.
- mtp reader: `mutdb.read_mtp_entry` (`mutdb.py:137`) is a **dispatcher** to
  `read_new_mtp_entry`; live callers at `mutdb.py:260` and `alchemy.py:1855`. Check
  whether `version='old'` is ever actually taken before touching it.

**Grab-bag modules (Phase 3 candidates):**
- `alchemy.py` (2547 lines): `mutate` router + `gen_hybrid_top`/`process_topol`
  (the topology-fill / gentop step) + a legacy inter-bonds recipe table +
  `_get_hybrid_residues`.
- `nsaa.py` (2379 lines): `mutate_nsaa` + the shared toolkit helpers + fragment
  machinery (`_prepare_nsaa_fragment`, `_find_methyl_cap`, …).

## Target shape
- **`covalent` = the single orchestrator + facade.** One `mutate`, one
  `build_hybrid_residue`, one `_build_target_residue`, one completion helper.
- **`nsaa.py` = the shared toolkit library** — keep the `_helpers`, drop the
  orchestration. (Not "legacy" — it stays the core; it just stops being an entry point.)
- **`nsaa.mutate_nsaa` → a thin shim** calling `covalent.mutate`. `pmx mutate_nsaa`
  CLI stays as an alias so nothing user-facing breaks.
- `alchemy.mutate` remains the classic API, delegating cleanly.

## Plan (risk-ordered; the test suite is the hard invariant)

### Phase 0 — golden-output guard  *(prerequisite for all phases)*
Snapshot the emitted artifacts for the matrix below, then require them
**byte-identical** through every refactor (plus `pytest -q` = 129 passing). The
matrix must span both **mutations and protonations** (protonation hybrids are how the
disulfide/NOS/pKa deprotonation legs are built — byte-identical there is
non-negotiable), plus each Stage-1 source and each special case:

| case | exercises |
|---|---|
| **Point mutations** | |
| ALA→VAL | small→β-branched, classic 20×20 |
| VAL→ALA | reverse (shrink → dummies) |
| ASP→LYS | charge-changing point mutation |
| LEU→PHE | **symmetric target (ring) — the MCS tie-break** |
| ILE→VAL | δ-atom / methyl symmetry (the historical ILE gotcha) |
| **Protonations** (same hybrid machinery, charge-changing) | |
| ASP↔ASPP (ASH) | carboxylate protonation |
| GLU↔GLUP (GLH) | carboxylate protonation |
| HSD/HSE→HSP | His tautomer → +1 (titratable neighbour coupling) |
| LYS↔LYN | ammonium ↔ neutral |
| CYS↔CYM | thiol ↔ thiolate (disulfide/NOS legs) |
| **Stage-1 sources / NSAA** | |
| LYS→ALY (structure template) | NSAA from structure, heavy-atom completion |
| LYS→ALY (from_rtp, no template) | template-free generative path |
| a `from_itp` ligand-fragment NSAA | ITP front end |
| **Terminal & bond** | |
| C-terminal amidation | terminal patch |
| N-terminal deletion | delete_terminal |
| a Cys–Cys disulfide | decouple_mutate / Stage-3 bond |
| a GAFF ligand covalent hybrid | ligand path + `gentop --bond` |

Capture per case: the hybrid `nsaa.rtp`/`nsaa.mtp` entry text, the `mutant.pdb`, and
the `gentop` output `.top`. Store as `tests/golden/` fixtures; add a test that diffs
them. Run each through **both** `pmx mutate` and `pmx mutate_nsaa` so the guard also
proves the two entry points agree before they're merged.
**Nothing else proceeds until this guard exists.**

### Phase 1 — kill the verified duplicates  *(low risk, high value)*
- Merge the two completion helpers into one shared function (natural home:
  `covalent/target.py`, since it's about building a target residue); `nsaa` imports it.
- Make `covalent._build_target_residue` the **single** target builder; have
  `nsaa.mutate_nsaa` call it instead of its inline load block. (Removes the fix-twice
  hazard permanently.)
- Confirm-then-merge (or leave with a note) the emission + mtp-reader items flagged
  above, gated on the golden guard showing identical output.

### Phase 2 — one orchestrator  ✅ DONE (2026-07-11)
- `nsaa.mutate_nsaa` is now a **thin adapter** over `covalent.mutate` (`--nsapdb` ->
  `structure`); every front end (rtp / structure / library / itp) goes through the single
  covalent orchestrator. The older 3-D `domapping` orchestration in `mutate_nsaa` is
  retired (~190 lines deleted); the dead `from_library` branch of `_build_target_residue`
  is removed. `pmx mutate_nsaa` CLI stays as a deprecated alias.
- **Correction to the plan's premise:** "Golden guard must stay byte-identical" was
  **infeasible** — the two orchestrators used different mappers (covalent graph-MCS vs
  nsaa 3-D `domapping`), which produce *different* (both valid) hybrids. Measured: LEU→PHE
  = identical mapping, only dummy *ordering* differs (ΔG-identical reindex); LYS→ALY =
  graph-MCS makes one fewer dummy H (28 vs 29 atoms — HN etc. mapped by identity, its
  claimed improvement). **User decision (confirmed): canonicalize on graph-MCS now; rebuild
  hybrids after.** The two structure/library goldens (`nsaa__lys2aly_struct`,
  `nsaa__leu2phe_lib`) were **re-captured** to the graph-MCS output; deterministic, 15
  golden cases, 0 drift, 130 tests. Any hybrid previously built via `pmx mutate_nsaa` must
  be rebuilt (topology reindexed; LYS→ALY-type cases change atom count).
- Remaining dead code (Phase 3): `nsaa.domapping` is now uncalled (kept — public-ish, part
  of the toolkit narrative; remove during module hygiene).

### Phase 3 — module hygiene  *(partly done 2026-07-11)*
- [x] **Dead `domapping` subsystem removed.** The 3-D-MCS "Atom mapping" section of
  `nsaa.py` (`domapping` + `_reformat_pdb`, `_set_backbone_carbonyl_double`,
  `_force_backbone_pairs`, ~190 lines) was uncalled after Phase 2 — deleted. Now-unused
  RDKit imports in `nsaa.py` reduced to a bare availability probe; stale `domapping`
  docstring refs (`covalent/target.py`, `covalent/mapping.py`) updated. `nsaa.py`
  −417 lines net this session. 18 golden cases, 0 drift, 130 tests.
- [~] **Keep `gen_hybrid_top`/`process_topol` in `alchemy.py`** (decision: moving them
  out leaves alchemy a near-empty shell). `alchemy.py` stays the router + emit/fill + FE
  analysis home.
- [ ] **DEFERRED: extract the shared toolkit from `nsaa.py` → `covalent/emit.py`**
  (`_make_rtp`/`_make_mtp`/`_merge_molecules`/dummy gen, the ~24 helpers `covalent`
  imports). High-churn, pure-structural (touches every import site); do as its own guarded
  step when there's appetite. `nsaa.py` currently = the shared toolkit + `mutate_nsaa`
  shim + the RDKit-free `sanitize_itp`/`prepare_nsaa_ff` ITP→FF path.

### Phase 4 — force-field data cleanup  *(DRAFT for review, 2026-07-11)*
**Why:** pdb2gmx globs *every* `*.rtp` in the FF dir, so a residue defined in two files
is a fatal `"Found rtp entries for 'X' in both …"`. The mutant FFs carry untracked runtime
leftovers that trigger exactly this (hit live this session: stale `K2AL.rtp` + runtime
`nsaa.rtp` both defining `K2AL`).

**Inventory (untracked leftovers, `git status` on each `*.ff`):**
| FF | untracked | per-hybrid `X2YY.*` | `nsaa.*` |
|---|---|---|---|
| `charmm36m-mut.ff` | 108 | **102** | 3 |
| `charmm36m-ecc-mut.ff` | 5 | 0 | 3 |
| `amber14sbmut.ff` | 2 | 0 | 2 |
| `amber14sb-q-mut.ff` | 4 | 1 | 0 |

- `X2YY.rtp/.arn` = **old per-hybrid emit convention** (one file per hybrid). Regenerable,
  and the collision source. → **delete.**
- `nsaa.{rtp,mtp,arn}` = **current** runtime accumulation (covalent engine writes here).
  Regenerable. → **delete from the tree + `.gitignore`** (don't commit runtime output).
- `mutres_term.{rtp,mtp,arn}` = **TRACKED / curated** (6 files, charmm36m-mut + -ecc). These
  are the pre-built terminal hybrid library — **KEEP** (verify they don't overlap nsaa.*).

**Proposed steps (each gated on your review — no deletion without confirming curated vs
regenerable):**
1. Settle the canonical emit convention = single per-FF `nsaa.{rtp,mtp}` (covalent already
   does this); retire the per-hybrid `X2YY.*` path entirely.
2. Delete the untracked `X2YY.*` and `nsaa.*` leftovers from all mutant FF dirs.
3. Add runtime-artifact globs to `.gitignore`: `*/nsaa.{rtp,mtp,arn}`, and the hybrid
   pattern (a `?2??`-stem matcher) so engine runs never re-dirty the tree.
4. Audit each shipped FF for duplicate residue defs across its committed `*.rtp`
   (`merged.rtp` vs `mutres.rtp` vs `mutres_term.rtp` vs `nsaa.rtp`) — the real pdb2gmx trap.
5. Regenerate the golden fixtures against the cleaned FF (should be 0-drift — the harness
   already strips these in its isolated copy, so the tree cleanup shouldn't move output).

### Phase 5 — FF-agnostic terminal framework  *(design, 2026-07-11)*
**Problem:** `_rtp_terminal` hardcodes CHARMM terminal patches (its docstring admits it), so
amber terminal deletion/mutation/amidation is broken (N valence-5; `KeyError 'OC'`).
**Key finding:** every patch it hardcodes is **already a named entry in the FF's own
`.n.tdb`/`.c.tdb`** — CHARMM ships `[ NH3+ ]`, `[ GLY-NH3+ ]`, `[ PRO-NH2+ ]`, `[ COO- ]`,
`[ CT2 ]` (the amide), and the `[ replace ]` blocks match `_rtp_terminal`'s `ren_atom`
calls line-for-line. The code re-transcribes FF data.
**Design:** replace the hardcoded patches with a generic **tdb applier**: parse the FF's
`.n/.c.tdb`, and apply the chosen variant's `replace`/`add`/`delete`/bonded ops to the RTP
entry; let RDKit place geometry from connectivity (so the tdb's geometry hints aren't even
needed). pmx's `term` arg maps to a variant name (`N`→`NH3+`, residue-selected `GLY-NH3+`/
`PRO-NH2+`; `C`→`COO-`; `C_amide`→`CT2`). Benefits: FF-agnostic; **consistent with pdb2gmx
by construction** (same tdb); the hardcoded PRO/GLY special cases **vanish** (they're tdb
entries). Catch: the pmx mutff **amber tdb files are EMPTY** → populate from the upstream
GROMACS amber ports, or support the alternate convention (separate `N<res>`/`C<res>` RTP
entries the amber rtp already carries). Fully generic form: `terminal_residue(res, ff, side,
variant)` = tdb-apply if tdb populated, else terminal-residue-variant lookup. Guard on both
CHARMM and amber. This folds `delete_terminal` / `amidate_cterminus` / terminal mutation
into one primitive (apply a terminus variant).

**Prototype validated (2026-07-11, `scratchpad/tdb_proto.py`).** A ~50-line tdb
parser + applier (ordered `replace`/`add`/`delete`/`impropers` ops; skip replace/add for
absent/existing atoms as pdb2gmx does; RDKit places geometry) diffed against the hardcoded
`_rtp_terminal` on CHARMM:
- N-term `NH3+` (ALA/LYS/VAL) and `GLY-NH3+` (GLY): **EXACT match.**
- C-term `COO-` (ALA/LYS): **EXACT match.**
- N-term `PRO-NH2+`: chemistry identical; only the two proton **names** differ (hardcoded
  `H1/H2` vs tdb `HN1/HN2`) — and the tdb form is what pdb2gmx emits, so the hand-coding is
  the odd one out.
- C-amide `CT2`: atoms/types/charges/bonds match; the tdb declares **2 impropers the
  hand-coding omits** (`C CA NT O`, `NT C HT2 HT1`) — a latent gap in `_rtp_terminal`.
Conclusion: the applier reproduces the hand-coded patches and is *more* faithful to pdb2gmx
where they differ.

**INTEGRATED into main (2026-07-11).** `nsaa._rtp_terminal(entr, term, ff=None,
resname=None)` now takes the patch from the FF's `.n/.c.tdb` when `ff`+`resname` are given
(parser `_parse_tdb`, applier `_apply_tdb_variant`, residue-prefix resolver
`_tdb_variant_name`: `GLY-NH3+`/`PRO-NH2+`/amber `LYS-NH3+` → generic `NH3+`/`COO-`/`CT2`),
falling back to the hardcoded CHARMM patches when no tdb. Threaded `ff`+`resname` through
`covalent.target._terminal_patch` (builds the terminal residue) and
`covalent.mutate.build_hybrid_residue` (entr1/entr2 impropers). Results:
- **CHARMM: byte-identical** for every terminal case except `CT2` amidation, which now
  gains the **2 impropers the hand-coding omitted** (`C CA NT O`, `NT C HT2 HT1`) — golden
  `term__cterm_amide` re-captured (correct improvement). 18 golden cases, 0 drift, 130 tests.
- **amber `-q`/`-q-ecc`: FIXED.** N/C-terminal deletion builds + pdb2gmx + gentop end-to-end
  on amber14sb-q / amber19sb-q (was the valence-5 crash). Per-residue variants resolve.
- **amber non-q (`amber14sb`/`19sb`): still the hardcoded fallback** (empty tdb) — the
  remaining gap. These ship explicit `N<res>`/`C<res>` RTP entries; the follow-up is the
  second resolver mechanism (use those entries when the tdb is empty).
Remaining follow-ups: (a) amber-non-q `N<res>`/`C<res>` resolver; (b) amber terminal golden
(parameterize the harness FF per case); (c) the hardcoded CHARMM fallback + its is_pro/is_gly
special cases are now superseded for tdb FFs — can be trimmed once amber-non-q is covered.

## Phase 6 — Boresch dummy-atom decoupling (`gentop --decouple-dummies`)  *(2026-07-12)*
**Why:** the Boresch audit (JCTC 2021, 17, 4403; paper `~/Downloads/dummy-atoms-...pdf`)
showed pmx's default "naive" dummy treatment keeps *every* bonded term of an
appearing/disappearing atom in BOTH end states, symmetrically. Some of those couple a
physical DOF to a dummy one at the junction and **break the factorization** that makes the
dummy free energy a state-independent constant (a systematic ΔG bias, common to classic and
generative). Flapping (#6–9) is ABSENT (pmx over-anchors). Correction is opt-in, off by
default, and must keep golden 0-drift when off.

**Design:** a new gentop-stage pass `alchemy._decouple_dummies(topol)`, gated behind
`gen_hybrid_top(..., decouple_dummies=)` → CLI `pmx gentop --decouple-dummies`. Runs after
the proline/terminal decouplings in `process_topol`. Off by default → golden byte-identical.

- **[x] Case #4 — Urey-Bradley terms (DONE, verified 2026-07-12).** CHARMM func-5 angles
  carry a UB (1-3 harmonic) term = a virtual bond between the two *outer* angle atoms.
  Boresch p.4406/4417: delete the UB "whenever the angle is made up from two physical and
  one dummy atom" (even the paper's *naive* baseline does this). A 2-phys/1-dummy angle has
  the dummy as an OUTER atom, so its UB 1-3 pair straddles the junction. Fix: zero the UB
  sub-term (S0,kUB) in the dummy state, KEEP the harmonic angle (theta0,ktheta) → no flapping
  introduced, no re-anchoring needed. **Aliasing gotcha:** `_find_angle_entries` stores the
  SAME list object in the A-slot `ang[4]` and B-slot `ang[5]` for a junction angle, so the
  pass *reassigns* the dummy-state slot with a fresh list (never mutates in place) — the
  physical state keeps its real UB, the dummy state loses it (correct asymmetry).
  Verified empirically (`scratchpad/verify_decouple.py`): on ala2val/asp2lys/leu2phe/
  ile2val/asp2aspp the ONLY lines that change are func-5 angle UB columns (`non-func5 = 0`),
  zeroed in exactly the dummy state; grompp-clean at both λ endpoints. Golden 0-drift (flag
  off), 130 tests. Note: 1-phys/2-dummy angle UBs (dummy central OR two outer dummies) are
  left per the paper's strict rule — they belong to the anchoring-count cleanup, not the
  crisp UB rule.
- **[x] Case #2/#3 — redundant terminal-junction dihedrals (DONE + validated 2026-07-13,
  `_dc_terminal_dihedrals`).** Scope (user-chosen): *strictly terminal, single-branch*
  junctions only — a dummy group bonded to the physical system through exactly ONE physical
  bridge atom, no ring back to physical, bridge carrying no second dummy group. Dual/triple/
  ring/shared-bridge are GUARDED OUT (ala2val/val2ala's two CB-methyls → 0 dihedral changes,
  UB only). Boresch p.4409 keeps *all angles* at a terminal junction; only dihedrals are
  redundant (ala2val: 23 junction dihedrals per methyl group; ile2val: 6 for the lone HV).
  Algorithm: pick one physical anchor atom (physical dihedral end-atom anchoring the most
  junction dihedrals; heavy-atom & low-id tie-break), KEEP its dihedrals unchanged, delete
  (zero, dummy state only, multiplicity preserved) every other junction dihedral into the
  group.
  **KEY EMPIRICAL LESSON — the paper's stiff re-anchor is WRONG for the coordinate-free
  gentop stage.** First cut followed Boresch literally (replace kept dihedral with a stiff
  mult-1 term, k=100 kcal/mol, Table 3 fn b). Two problems, both caught by validation:
  (a) GROMACS forbids perturbing a proper-dihedral *multiplicity* between states A/B
  (`multiplicity can not be perturbed 1!=3`); (b) even after working around (a) by appending
  a separate mult-1 entry, an anchor-fc sweep (`scratchpad/`) showed the stiff anchor at a
  FIXED phase (can't match the absent geometry) COUPLES the physical bridge to the dummy and
  displaces physical HEAVY atoms up to 244 pm (leu2phe) / 46 pm (ile2val) — worse than the
  redundancy removed. Deletions ALONE: <1–5 pm. **Fix: keep the anchor atom's REAL dihedral
  untouched** (its periodic minima already match the built geometry; the dummy's torsion
  profile then equals the physical atom's), no artificial stiff term. Validated: physical
  heavy-atom RMSD (naive vs corrected, vacuum min at dummy endpoint) 0.35–3.0 pm (ile2val
  0.75, leu2phe 3.0, ala2pro 0.35, asp2lys 1.6); gas-phase flapping MD
  (`scratchpad/flapping_check.py`, 100 ps) — pruned+kept-real group localizes as well as the
  over-anchored naive (unimodal, antipodal excursions ≤ naive) for single-H, methyl, and
  PHE-ring groups; grompp-clean both λ endpoints; golden 0-drift (flag off); 130 tests.
  NOT done: dual/triple/ring junctions (guarded out; need the fuller Boresch §2.3 pass).
- **[x] Factorization spot check for case #4 (DONE 2026-07-12, `scratchpad/factorization_check.py`).**
  Structural form of the paper's check #1: at the endpoint where the group is a dummy,
  vacuum-minimize with the naive vs corrected topology and measure the PHYSICAL-atom
  displacement (the clean factorization-break signal; total Epot is contaminated by the
  nonbonded-inert dummy group's own bonded relaxation, so it is NOT reported as the bias).
  Marginal effect of the UB removal on physical geometry (tight-converged):
  | case (endpoint) | phys-atom RMSD | max displacement |
  |---|---|---|
  | ILE→VAL, δ-methyl dummy (B) | 0.66 pm | 1.7 pm |
  | ALA→PRO, ring appears (A) | 0.49 pm | 2.2 pm |
  | ASP→LYS, NZ amine + carboxylate (A) | ~2.0 pm | ~7.5 pm |
  Conclusion: **case #4 in isolation is a small, real correction** — pmx over-anchors, and
  the retained harmonic angle stays as the anchor, so removing the UB shifts the physical
  minimum only sub-few-pm. Largest for stiff heteroatom junctions (Asp carboxylate kUB up to
  1.76e5, Lys amine); negligible for a plain methyl/ring. The theoretical value is restoring
  *exact* factorization (a systematic bias, small per-junction, additive over many junctions);
  the bulk of the remaining factorization break is in the un-pruned redundant dihedrals (#2/#3).

## Guardrails
- **Behavior-preserving only.** Golden fixtures + 129 tests gate every commit; if any
  byte changes, stop and explain before proceeding.
- **Docs/pmx only** — the `GitHub/pmx` checkout is off-limits.
- Update `residuetypes.dat` registration consistently (both orchestrators touch it).
- After Phase 1, migrate the showcase build scripts from `pmx mutate_nsaa` to
  `pmx mutate` so tooling sits on the canonical entry.

## Recommendation
Do **Phase 0 + Phase 1** first — it removes the concrete "fix-it-twice" hazard with
near-zero risk. Treat Phases 2–3 as separate, schedulable follow-ups.

---
## Phase 0 status: BUILT (2026-07)
Golden fixtures + guard live in `tests/golden/` (`generate_golden.py`, `test_golden.py`,
shared `_peptide.pdb`). **130 tests pass** (129 + the golden guard). 13 cases captured,
all regenerate **byte-identical** (0 drift):
- 11 full (mutant.pdb + gentop `hybrid.top`): 5 point mutations, 5 protonations
  (ASPP/GLUP/HSP/LYN/CYM), C-terminal amidation.
- 2 partial (mutant.pdb only): the NSAA from-structure and from-rtp cases — these
  **surfaced a real bug** the guard is meant to catch: `mutate_nsaa` raises a
  non-exception (`TypeError: exceptions must derive from BaseException`) and the
  from_rtp generative path fails downstream pdb2gmx. Fix those, then their `hybrid.top`
  completes the matrix. Generative cases are isolated in a temp FF copy (verified: the
  tracked FF is not polluted).
Notes: fixtures generated against the **canonical** `src/pmx/data/mutff` FF (the
showcase-bundle FF is mutated and gave false gentop failures — do not use it here).

## Phase 1 progress (2026-07)
- [x] **Completion helper unified.** `nsaa._complete_target_from_rtp` +
  `covalent._complete_target_from_structure` (near-identical Kabsch-overlay) → a single
  `covalent.target.complete_target_from_structure`, imported by both. Removed the
  duplicated ~30-line body; the "structure = heavy-atom geometry, RTP = identity/types/
  H" rule now has one home. **Golden guard 0 drift, 130 tests pass.**
  This proves the safe-refactor loop (edit → golden --check → suite).
- [x] **Target *builder* — structure front end unified.** `nsaa.mutate_nsaa`'s
  pure-structure branch (`nonstandardPDB`, no `itp`) now delegates to the single
  `covalent._build_target_residue(..., structure=...)` — the same builder the covalent
  orchestrator uses. This is the real fix-twice root (the heavy-atom bug lived here in
  both paths). **Golden 0 drift, 130 tests pass**, and it is genuinely guarded: the
  `nsaa__lys2aly_struct` fixture now carries a full `hybrid.top`, so any drift in the
  emitted dual-state topology fails the tripwire. (Minor, guard-confirmed harmless deltas
  vs the old inline block: the builder `set_resname(target)`, `attach_bonds_from_rtp`, and
  assigns RTP params unconditionally rather than gated on `resname in rtp` — the last turns
  a "structure target absent from RTP + no itp" misuse from a warn-and-limp into a clear
  early error, matching covalent.)
- [x] **`from_library` (build_chain) branch unified.** The single builder gained a
  `from_library=` front end that keeps the internal-library geometry (NOT `from_rtp` — the
  two give the same atom set but different atom *order* for e.g. TRP, so `from_rtp` would
  drift). `nsaa.mutate_nsaa` now routes **every non-ITP** front end (structure *and*
  from_library) through `covalent._build_target_residue`; only the capped-ITP branch stays
  inline. Guarded by a new golden case **`nsaa__leu2phe_lib`** (LEU5→PHE via build_chain +
  the nsaa `domapping` engine; full `hybrid.top`) — deterministic, **0 drift, 130 tests**.
  Net: `_build_target_residue` is now the single Stage-1 target builder for the RTP,
  structure, and library front ends, shared by both orchestrators.
- [x] **`itp` (capped-fragment) branch unified — `pmx mutate_nsaa --itp` is now a thin
  adapter over `pmx mutate --itp`.** Decisions taken (user-confirmed): (a) a structure is
  now **optional** with `--itp` — geometry embeds from the ITP bond graph via
  `TargetResidue.from_itp`, `--structure`/`--nsapdb` only overrides; (b) the covalent
  graph-MCS engine is the **single** ITP implementation — `mutate_nsaa`'s bespoke
  structure-required + `domapping` ITP path is retired; (c) `mutate_nsaa(itp=...)` now
  early-returns by delegating to `covalent.mutate(structure=nonstandardPDB, itp=...)`, so
  the hybrid accumulates into `nsaa.{rtp,mtp}` (auto-discovered — no `--extra_mtp`).
  Verified `mutate_nsaa(itp, nsapdb=None)` is **byte-identical** to `covalent.mutate(itp)`.
  New golden **`nsaa__leu2sfr_itp`** (LEU5→SFR via a CH3-capped fragment `sfr_frag.itp`, no
  structure, embed path) exercises the delegation end-to-end; deterministic, full
  `hybrid.top`. CLI help in `scripts/mutate_nsaa.py` updated (structure optional; ITP path
  reports the no-`--extra_mtp` behaviour). **15 golden cases, 0 drift, 130 tests.**
  Net: `_build_target_residue` + `covalent.mutate` are the single Stage-1 builder /
  orchestrator for ALL front ends (rtp, structure, library, itp); `mutate_nsaa` keeps its
  own `domapping` engine only for the structure/library paths.
