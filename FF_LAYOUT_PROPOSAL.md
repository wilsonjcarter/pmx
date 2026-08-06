# Force-field layout proposal — separating introduced residues and parameters

Status: **proposal, nothing implemented.** Written 2026-08-06 after the NOS/CN0 loss.

---

## 1. Why the current layout keeps failing

New pmx builds hybrids on the fly *from RTP entries*, so an RTP entry is no longer
documentation — it is the input to the engine. But every residue and parameter we have
introduced currently lives **inside a vendored file**:

| what we introduced | where it lives today | tracked? |
|---|---|---|
| `GACE` cap (GSH reference) | `amber*-q-ecc-mut.ff/aminoacids.rtp`, `charmm36m-ecc-mut.ff/merged.rtp` | tracked |
| `ALY` acetyl-lysine | `amber14sb-q-ecc-mut.ff/aminoacids.rtp:544` (**indented**, inserted mid-file) | tracked |
| `NME` cap | `charmm36m-ecc-mut.ff/merged.rtp:17` (**indented**) | tracked |
| sulfenic/NOS bonded terms | `ffbonded.itp` tail | tracked |
| protein↔CGenFF link terms | `cgenff_link.itp` | **untracked** |
| `CN0`/`KN0`/`KS0` endpoints | *(were in `merged.rtp`)* | **lost** |
| generated hybrids | `nsaa.{rtp,mtp,arn}` | **gitignored by design** |

Three distinct failure modes follow from this, and we have now hit all three:

1. **Revert/regeneration clobbers.** `merged.rtp` and `ffbonded.itp` are tracked *vendored*
   files. A `git checkout` or a regeneration pass silently discards hand-added entries
   because they are indistinguishable from vendor content. This is exactly how CN0/KN0/KS0
   and the sulfenic bonded block disappeared between 2026-07-13 and 2026-07-15.
2. **Untracked survives but is invisible.** `cgenff_link.itp` outlived the revert precisely
   because git does not manage it — but it is also not reviewed, not diffed, and would die
   to `git clean -fdx`. It was present but not `#include`d, so it did nothing.
3. **No way to tell ours from theirs.** There is no diff, no manifest, and no test that
   answers "which residues did we add, and are they all still here?"

The 2026-08-05 recovery only worked because the Desktop bundle
`~/Desktop/pmx_covalent_showcases/ff/` happened to be a snapshot from before the loss, and
because built topologies carry explicit per-interaction parameters. Neither is a backup
strategy.

---

## 2. The constraint that decides the design

**pdb2gmx globs `*.rtp` in the FF directory. pmx does not.**

pmx resolves an RTP residue through a single hardcoded path, repeated in four places:

| call site | what it does |
|---|---|
| `src/pmx/covalent/target.py:148` | `TargetResidue.from_rtp` — the generative engine's entry point |
| `src/pmx/covalent/mutate.py:91` | hybrid build |
| `src/pmx/covalent/mapping.py:130` | graph mapping |
| `src/pmx/nsaa.py:1466`, `:1692` | NSAA emit path |

all of them:

```python
rtp_file = os.path.join(ffpath, "merged.rtp")
if not os.path.isfile(rtp_file):
    rtp_file = os.path.join(ffpath, "aminoacids.rtp")
```

So **splitting residues into topic files works for pdb2gmx immediately and silently breaks
pmx.** Any split must be preceded by unifying these four call sites behind one loader that
globs. This is the whole reason the split is a code change and not a file move.

There is already a precedent for the glob-and-merge pattern: `alchemy.py:2340` auto-discovers
`nsaa.mtp` and folds its names into `extra_name_map`. The RTP loader should work the same way.

Two hard rules the split must respect, both learned the hard way yesterday:

* **Move, never copy.** GROMACS rejects a residue defined in two `.rtp` files, and grompp
  rejects "a second block of parameters for the same atom types" — appending the bundle's
  bonded block on top of ECC's existing copy broke three previously-working legs. Entries
  must be *removed* from the vendored file in the same commit that adds them to the new one.
* **Precedence must be explicit.** If a topic file and a vendored file both define a residue,
  the loader needs a defined winner, not dict-update order.

---

## 3. Proposed layout

One file per *chemistry*, per force field, sitting alongside the vendored files:

```
<ff>.ff/
  merged.rtp / aminoacids.rtp     vendored, NEVER hand-edited again
  ffbonded.itp / ffnonbonded.itp  vendored, NEVER hand-edited again
  pmx_ptm.rtp                     phospho + methyl + acyl  (SEP TPO PTR / MLY M3L / ALY ...)
  pmx_redox.rtp                   sulfenic + NOS + disulfide endpoints (CSO CN0 KN0 KS0 CYS2 CYM)
  pmx_caps.rtp                    GACE, NME, and any other cap we introduced
  pmx_extra_bonded.itp            every bonded term we added, #include'd from pmx_extra.itp
  pmx_extra_nonbonded.itp         every atomtype/nonbonded term we added
  pmx_extra.itp                   the single #include hook (pulls the two above + cgenff_link)
  pmx_extra.manifest              machine-readable list of what we introduced (see §5)
  nsaa.{rtp,mtp,arn}              generated, stays gitignored
```

Naming: a `pmx_` prefix makes provenance obvious in a directory listing and makes
"did we add this?" a glob rather than an archaeology session.

**Only one new `#include` in a vendored file** — `forcefield.itp` gains `#include
"pmx_extra.itp"`, and everything else hangs off that. That single line is the only vendored
edit in the design, so a regeneration has exactly one thing to re-apply, and a missing
`pmx_extra.itp` is a loud failure rather than a silent one.

### Which residues move, by force field

Grounded in the inventory (`python3` scan of all 16 `*.ff` dirs):

* **CHARMM** (`charmm36m-mut`, `charmm36m-ecc-mut`) — most PTMs genuinely ship with CHARMM:
  SEP/TPO/PTR, MLY/M3L and CSO are all real vendor entries in `merged.rtp`. **Leave them.**
  Only the things we introduced move: `NME` (indented at `merged.rtp:17`), `GACE`, `TYM`,
  and the CN0/KN0/KS0 we are about to create → `pmx_redox.rtp` / `pmx_caps.rtp`.
  Non-trivial call: CSO stays vendored but its *bonded terms* are ours, so the residue and
  its parameters end up in different files. That is correct — provenance follows origin,
  not topic — but it must be documented or it looks like an oversight.
* **Amber** (`amber14sb*`, `amber19sb*`, `amber99sb*`) — the opposite: the vendored files
  have **no** phospho, methyl or acyl entries at all. Everything of that kind is ours.
  `GACE` and `ALY` currently sit inside `aminoacids.rtp` and move out wholesale into
  `pmx_ptm.rtp` / `pmx_caps.rtp`. This is where the split buys the most.
* `charmm22star-mut`, `amber99sbdisp`, `des-amber` — nothing introduced; no new files.
  Absence of `pmx_*` files is itself informative.

---

## 4. Backup, before anything moves

1. **Snapshot the whole tree first**, to the crucial drive (main drive is at 86–97%):
   `/Volumes/crucialP310/mutff_snapshots/mutff_2026-08-06_pre-split/`. Cheap, and it is the
   thing we did not have on 2026-07-15.
2. **Fold in the Desktop bundle.** `~/Desktop/pmx_covalent_showcases/ff/` currently holds
   content the core tree has lost — the 212-entry `nsaa.rtp` (incl. `C2CN`/`L2KN`/`L2KS`)
   and the sulfenic block. Snapshot it too and treat it as the recovery source for CN0/KN0.
3. **Stop relying on `nsaa.*` being disposable.** It is gitignored on the theory that it is
   regenerable, which is true only while the *inputs* survive. Once the endpoints live in a
   tracked `pmx_redox.rtp`, that theory actually holds.

---

## 5. Making the loss detectable

The split is only half the fix. The other half is that **nothing today notices when an entry
disappears.** Proposed, cheapest first:

* **`pmx_extra.manifest`** — one line per introduced residue/parameter block with the file it
  belongs in. A test walks every `*.ff` dir, reads the manifest, and asserts each entry
  resolves. Turns a silent revert into a failing test.
* **Extend the golden guard** (`tests/golden`) with one build per chemistry — a phospho, an
  acyl, a sulfenic, and a NOS bond-forming leg. Any of yesterday's failures would have been
  caught the moment it was introduced.
* **A `pmx ff doctor` subcommand** — for a given ff, report which `pmx_*` files are present,
  whether `forcefield.itp` includes `pmx_extra.itp`, and which manifest entries fail to
  resolve. This is what I would have run yesterday instead of an hour of archaeology.

---

## 6. Suggested order

1. Snapshot (§4) — do this regardless of what else gets decided.
2. Unify the four RTP call sites behind one glob-aware loader, with precedence rules and
   tests. **No files move in this step**; behaviour must be identical afterwards.
3. Reconstruct `CN0`/`KN0`/`KS0` into `pmx_redox.rtp` for `charmm36m-ecc-mut`, derived from
   the bundle's `C2CN`/`L2KN`/`L2KS` hybrid entries. Validate by rebuilding a NOS
   `nos_formation` leg and comparing against the Jul-13 topology.
4. Move CHARMM's introduced entries out (small, low risk — mostly `NME`/`GACE`/`TYM`).
5. Move Amber's PTM/cap entries out (larger, and the real prize).
6. Separate the bonded/nonbonded additions into `pmx_extra_*.itp` — **move, not copy**.
7. Manifest + golden-guard coverage + `pmx ff doctor`.

Steps 1–3 unblock the NOS showcase on their own; 4–7 are the durable fix.

---

## 7. Open questions

1. **Per-ff or shared topic files?** The proposal is per-ff (parameters are ff-specific, so
   sharing means a lookup layer). A shared `pmx_ptm.rtp` across the six Amber variants would
   cut duplication a lot — but only if their atom types genuinely agree. Worth checking
   before committing to per-ff duplication.
2. **How far does "audit" go?** The 16 `*.ff` dirs also differ in ways unrelated to this:
   `charmm36m-ecc-mut` is missing `mutres_term.*` and `residuetypes.dat` that
   `charmm36m-mut` has, and has extra `aminoacids.r2tp` / `extract_canonical_aas.py`. Is
   normalising that in scope, or strictly the introduced-residue split?
3. ~~Do the ECC variants need their own PTM parameters?~~ **ANSWERED — and it is a live bug.
   See §8.**
4. **What is the vendored baseline?** To ever diff "ours vs theirs" mechanically we need a
   pristine copy of each upstream ff. Do we have one, or does it need fetching?

---

## 8. Found while scoping this: the ECC PTM residues are NOT scaled

Net charges via `pmx.forcefield.RTPParser` on `merged.rtp`:

| residue | `charmm36m-mut` | `charmm36m-ecc-mut` | verdict |
|---|---|---|---|
| LYS | +1.0 | **+0.75** | scaled ✅ |
| ARG | +1.0 | **+0.75** | scaled ✅ |
| ASP | −1.0 | **−0.75** | scaled ✅ |
| GLU | −1.0 | **−0.75** | scaled ✅ |
| HSP | +1.0 | **+0.75** | scaled ✅ |
| CYM | −1.0 | **−0.75** | scaled ✅ |
| **SEP** | −1.0 | **−1.0** | **NOT scaled ❌** |
| **TPO** | −2.0 | **−2.0** | **NOT scaled ❌** |
| **PTR** | −1.0 | **−1.0** | **NOT scaled ❌** |
| **MLY** | +1.0 | **+1.0** | **NOT scaled ❌** |
| **M3L** | +1.0 | **+1.0** | **NOT scaled ❌** |

The ECC pass covered the **canonical** charged residues and stopped there. Every
phosphorylated and methylated residue in `charmm36m-ecc-mut` still carries full formal
charge, so any ECC phospho/methyl calculation mixes scaled protein charges with an unscaled
modification — `TPO` at −2.0 being the worst case. Neutral PTMs (`CSO`, `ALY`, `LSN`,
`GACE`, `NME`) are unaffected, since there is nothing to scale.

This is precisely the class of error the split is meant to prevent: the entries look
authoritative, sit in a vendored file, and nothing checks them. It also changes the order of
work — **scaling must be decided before the PTM entries are moved**, or the split will
enshrine the wrong charges in a file that reads as curated.

Which ECC charges are *correct* for a phosphate is a real question, not a mechanical ×0.75:
the scaling factor depends on the ionic group, and a −2 phosphate is not obviously
`2 × 0.75`. This needs the ECC reference, not an analogy.

---

## 9. Upstream baseline found — and what it changes

`/Users/carter/Documents/GitHub/pmx` **is** upstream `wilsonjcarter/pmx` @ `develop`
(5061936b), so open question 4 is answered: we have a pristine baseline and "ours vs theirs"
is now a mechanical diff. (Read-only — that checkout stays off-limits for writes and imports.)

### Ours (present here, absent upstream)

| ff | files |
|---|---|
| `amber14/19sb-q-mut`, `amber19sbmut`, `amber99sb-star-ildnp-q-mut` | `disulfide.{rtp,mtp,hdb}` |
| `charmm36m-mut`, `charmm36m-ecc-mut` | `cgenff_link.itp`, `nsaa.{rtp,mtp,arn}` |
| `charmm36m-mut` | `opc.itp`, `residuetypes.dat` |
| `amber14sb-q-mut` | `Archive.zip` ← junk, delete |
| `charmm36m-ecc-mut` | `TYM` in `merged.rtp` |

Content-wise `charmm36m-ecc-mut/merged.rtp` is upstream **+ TYM and nothing else** — which
confirms the CN0/KN0/KS0 diagnosis: they were local additions to a tracked vendored file, and
the file is now back at its upstream state. Nothing was "lost from upstream"; our edit was
reverted. `NME` turns out to be genuine upstream content, not ours — I had inferred otherwise
from its indentation, which was wrong.

### Lost from here (present upstream, absent now)

| ff | files | note |
|---|---|---|
| `amber14sb-q-mut` | `S2ER.rtp S2P1.rtp S2P2.rtp S2SE.rtp S2SP.rtp` | **phospho** hybrids (SER→SEP etc.) |
| `charmm36m-mut` | `L2YS.{arn,rtp} S2LY.{arn,rtp}` | LEU→LYS, SER→LYS hybrids |
| `charmm36m-ecc-mut` | `mutres_term.{arn,mtp,rtp}` | terminal patches |

**A second silent-loss vector, same class as the first.** `git check-ignore` confirms
`.gitignore:51-53` (`?2??.rtp` / `?2??.mtp` / `?2??.arn`) matches `S2P1.rtp`, `L2YS.rtp` and
friends. The pattern was written to keep *generated* single-hybrid files out of the repo, but
upstream **ships some of them as curated content** — so locally they are invisible,
unrestorable by git, and deleted without trace. The gitignore needs to distinguish generated
from curated, not pattern-match on filename shape.

### The linchpin bug, now demonstrable

`disulfide.rtp` already **is** the proposed layout — a tracked, per-topic hybrid file, shipped
for four Amber force fields. But it holds hybrids (`[ C2CD ] ; Cys -> Cym`), while
`from_rtp` needs the *endpoint* residue. Checking `CYS2` visibility:

| force field | `CYS2` in main rtp | has `disulfide.rtp` |
|---|---|---|
| `charmm36m-ecc-mut` | ✅ | — |
| `charmm36m-mut` | ✅ | — |
| `amber14sb-q-mut` | ❌ | ✅ |
| `amber19sb-q-mut` | ❌ | ✅ |
| `amber19sbmut` | ❌ | ✅ |
| `amber99sb-star-ildnp-q-mut` | ❌ | ✅ |

So **the disulfide recipe works on CHARMM and would fail on all four Amber force fields**,
with the identical `KeyError: Residue 'CYS2' not found in .../aminoacids.rtp` that CN0 threw.
It has not been hit only because the disulfide showcase runs on CHARMM. This is a latent bug
today, independent of any reorganisation, and it makes the glob-aware loader (§2, step 2) the
linchpin rather than a nicety.

### Revised naming

Given `disulfide.{rtp,mtp,hdb}` already exists and is tracked, drop the invented `pmx_`
prefix and follow the established convention: `phospho.rtp`, `methyl.rtp`, `acyl.rtp`,
`sulfenic.rtp`, `nos.rtp`, `caps.rtp`. Provenance then comes from the manifest (§5) and the
upstream diff, which is more reliable than a filename prefix anyway.

### Parsing note

Three separate times while auditing, a naive `^\s*\[\s*NAME\s*\]\s*$` regex missed real
entries: `[ NME ]` (leading space), `[ ALY ]` (indented, mid-file), `[ C2CD ] ; Cys -> Cym`
(trailing comment). Any manifest/doctor tooling must parse via `RTPParser`, not grep — and
the fact that these files resist mechanical auditing is itself part of the problem being
fixed.

---

## 10. Implemented 2026-08-06

Decisions taken: **per-force-field** topic files (not shared); **ECC PTM scaling deferred**
to its own workstream; disulfide needs only the *endpoint* entry — hybrids stay in the legacy
`mutres.rtp`, since pmx generates `nsaa` RTPs on the fly now.

### Backup

`/Volumes/crucialP310/mutff_snapshots/2026-08-06_pre-split/` — 446 MB, three trees:
`core_mutff` (working), `desktop_bundle_ff` (only surviving copy of the lost CN0/KN0
hybrids + sulfenic block), `upstream_mutff_5061936b` (pristine baseline). `README.txt`
records what each is.

### Code: one glob-aware RTP loader

`ffparser.load_ff_rtp(ff)` parses the base RTP and merges every non-hybrid `*.rtp` in the ff
directory, exposing `.sources` and `.conflicts`. Base wins a clash; clashes are recorded, not
raised, so an unsplit force field behaves exactly as before. `_is_hybrid_rtp` excludes the
hybrid libraries (`mutres*.rtp`, `nsaa.rtp`, `disulfide.rtp`) and the `?2??.rtp` single-hybrid
pattern.

All five hardcoded call sites now route through it: `covalent/target.py` (`from_rtp`, still
`lru_cache`d — now keyed on the ff rather than a file path), `covalent/mutate.py`,
`covalent/mapping.py`, `nsaa.py` ×2. Lookup failures name the directory and the files
actually searched instead of a single path.

Verified as a no-op on the current tree: every force field reports the same entry count
before and after, no topic files existed yet to merge.

### Amber disulfides were broken; now fixed

`_RESNAME_SYNONYMS` mapped `CYX -> [CYX, CYS2]` but not the reverse, so a recipe naming the
`CYS2` endpoint failed on every Amber force field. Added `"CYS2": ["CYS2", "CYX"]`. `CYS2`
now resolves to `CYX` (10 atoms, matching CHARMM's `CYS2`) on `amber14sb-q-mut`,
`amber19sb-q-mut`, `amber19sbmut` and `amber99sb-star-ildnp-q-mut`.

### Eight lost residues restored as topic files

Diffing every base RTP against upstream: **only `charmm36m-mut` had lost content** — 1047
entries upstream, 1039 here. Every other force field had only gained entries. Restored
*verbatim* (raw source text, comments and charge formatting preserved — not re-serialised
through `RTPParser`):

* `charmm36m-mut.ff/acyl.rtp` — `ALY BLY FLY HIB KCR PLY SLY` (the acyl-lysine set)
* `charmm36m-mut.ff/caps.rtp` — `NME`

The force field is back to 1047 entries, zero conflicts. `NME`'s absence is exactly what made
the GACE reference arm fail pdb2gmx with `Residue type 'NME' not found` during the 2026-08-05
NOS work — that failure is now explained and fixed.

### Golden guard

`nsaa__lys2aly_struct` and `nsaa__lys2aly_rtp` had been failing because `ALY` was gone.
Both pass now: **16 ok, 0 ERROR**. 129 unit tests pass.

Two pre-existing failures remain and are untouched: `term__cterm_del` / `term__nterm_del`
(`no-top`, `hyb.top` never produced). Present in every run including the HEAD baseline.

### The gitignore question — resolved as a non-issue

`?2??.rtp` (`S2P1.rtp`, `L2YS.rtp`, ...) *are* generated hybrids, so the ignore rule is
correct and no restore is needed: pmx regenerates them on the fly. §9's "second silent-loss
vector" was wrong on that point.

### Open, deliberately not done

`charmm36m-ecc-mut` is missing `mutres_term.{arn,mtp,rtp}` that upstream ships — an unstaged
tracked deletion, restorable with `git checkout --`. Left alone: the ECC tree shows signs of a
deliberate regeneration (`aminoacids.r2tp`, `extract_canonical_aas.py`), so the deletion may
have been intentional. **Needs a decision.**

Next: reconstruct `CN0`/`KN0`/`KS0` into `charmm36m-ecc-mut.ff/nos.rtp` from the bundle's
`C2CN`/`L2KN`/`L2KS` hybrids — the loader now makes such a file visible to the engine.

---

## 11. CN0 / KN0 / KS0 reconstructed and validated (2026-08-06)

**All 7 NOS legs build again** — the full recipe has not worked since 2026-07-13.

### Reconstruction

The endpoints were not re-derived, they were read back out of runs that worked:

| residue | source | atoms | qnet |
|---|---|---|---|
| `CN0` | state-B columns of `C2CN` in `CNOS/fe_calcs/1e8i/nos_formation/pmxtop.top` | 11 | −0.154 |
| `KN0` | state-B columns of `L2KN`, same topology | 20 | +0.154 |
| `KS0` | state-B columns of `[ L2KS ]` in the bundle `nsaa.mtp` | 19 | +0.286 |

Bonds/impropers/cmap come from the same hybrids, restricted to the state-B branch of each
split atom and renamed (`SG2 -> SG`, `CB2 -> CB`, `NZ2 -> NZ`, ...). Written to
`charmm36m-mut.ff/nos.rtp`.

**CN0 + KN0 = 0.000 exactly.** The two halves were extracted independently; that the bridged
pair comes out neutral, with ±0.154 transferred across the S–O⋯N bridge, is the check that
the extraction is right.

The bundle's `nsaa.rtp` turned out to be **damaged** — lines 155–371 are orphaned
`[ atoms ]`/`[ bonds ]` bodies with no residue header, which is why `RTPParser` reported
`C2CN` with 92 atoms (it had swallowed four following entries). Reading the state-B columns
out of a *topology* rather than that file avoided the corruption entirely.

### Validation against the July reference

Rebuilt `1vhf` and compared the emitted hybrid to `CNOS/fe_calcs/1vhf/nos_formation/pmxtop.top`
(same system, 2026-07-13):

* `L2KN` — **26/26 atoms identical** (types and charges, both states).
* `C2CN` — **14/17 identical**. The 3 differences are only the *state-A* (CSO) `CB`/`HB`
  atom **types**: `CG321/HGA2` now vs `CT2/HA2` in July. All charges identical; the state-B
  half, which *is* CN0, matches on every atom.

That difference is deliberate. July's `merged.rtp` had CSO locally retyped to protein types,
and that edit was reverted along with everything else. Rather than re-patch a vendored
residue (which the loader would not let a topic file override anyway), the CGenFF-typed
bonded terms are supplied in `nos_bonded.itp`.

### Three more losses found and fixed along the way

1. **`charmm36m-mut/ffbonded.itp`** had lost the whole NOS/sulfenic bonded block *and* never
   included `cgenff_link.itp`. Restored into `nos_bonded.itp`, `#include`d once from the end
   of `ffbonded.itp`.
2. **`NME`'s hydrogen-build entry** was missing from `merged.hdb`, so pdb2gmx found `[ NME ]`
   but could not build its hydrogens → `4 missing atoms in molecule Protein`. Restored.
3. **`charmm36m-mut.ff/residuetypes.dat` was a 12-byte stub containing only `6FL Protein`.**
   A per-force-field `residuetypes.dat` *replaces* the root one rather than merging, so this
   silently stripped `ACE`/`NME`/`CN0`/`KN0` of their `Protein` classification for this force
   field alone. `6FL` was **already** in the root file, so the stub was pure harm. Moved aside
   to `residuetypes.dat.bak_6FL_only`; added the missing `KS0 Protein` to the root file.

### A constraint the topic-file design has to respect

`NME` **does not work from a topic file.** pdb2gmx globs `*.rtp` and `*.hdb` correctly (the
log confirms it opens `acyl.rtp`, `caps.rtp`, `nos.rtp`, `caps.hdb`), and it resolves the
entry — but a terminal cap resolved out of a secondary RTP fails at atom-sorting with
`Atom H in residue NME 5 was not found in rtp entry NME with 6 atoms`. Moving the identical
entry into `merged.rtp` fixes it; nothing else changed. Root cause not established.

This gives a cleaner rule than the original proposal:

> **Upstream content belongs in the vendored file; topic files are for what we add beyond
> upstream.** The 8 lost `charmm36m-mut` entries were upstream content, so restoring them to
> `merged.rtp` re-establishes parity — and a future revert *to upstream* then keeps them
> rather than dropping them. `CN0`/`KN0`/`KS0` exist in no upstream, so they are exactly what
> a topic file is for, and `nos.rtp` is validated working.

Current split for `charmm36m-mut`: `merged.rtp` +NME (upstream parity), `acyl.rtp` (the 7
acyl-lysines — upstream content, but working from a topic file and matching the intended
per-chemistry layout), `nos.rtp` + `nos_bonded.itp` (ours).

`merged.rtp` must **not** be replaced with upstream wholesale: 73 lines exist only here,
explicit `TPO` bonds/angles from the RTPParser angles fix. The merge has to stay surgical.

### Regression

129 unit tests pass. Golden guard: **16 ok**, plus the same 2 pre-existing `no-top`
(`term__cterm_del`, `term__nterm_del`) present in every run including the HEAD baseline.

### Not done

`nos.rtp` was **not** copied into `charmm36m-ecc-mut` — its charges are the unscaled
`charmm36m-mut` values, and dropping them into the ECC force field would repeat the mistake
documented in §8. That waits on the ECC scaling decision.
