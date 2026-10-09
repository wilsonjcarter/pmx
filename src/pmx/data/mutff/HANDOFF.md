# Clean mutff — handoff (continue here)

Goal: `src/pmx/data/mutff_clean/` = **pristine OG base** (upstream/deGrootLab `develop`) + pmx
additions as **topic files only**. Base files never edited (except NME, see below). See
`CLEAN_FF_PLAN.md` for the design. OG `mutff/` is untouched; nothing deleted; all work is new.

## Done
- **charmm36m-mut.ff**: OG base + NME wired into base `merged.rtp` (after ACE) + hdb. Topic files:
  `acyl.rtp`, `nos.rtp`+`nos_bonded.itp`+`cgenff_link.itp`, `mutres_term.*`, `ptm_phospho.rtp`
  (SEP/TPO/PTR/YP1/YP2), `ptm_methyl.rtp` (MLY/**MLYS**/M3L), `ptm_misc.rtp` (CSO/CSS/NEP/ASPM),
  `glycation.rtp` (KCX), `ptm_extra.rtp` (CIR/HYP/MLZ/SNC). Atom-type coverage clean.
- **amber14sbmut.ff + amber99sb-star-ildn-mut.ff**: OG base + topic files:
  `ptm_phospho.rtp` (**SP1/SP2/TP1/TP2/YP1/YP2** — the engine's target names, NOT SEP/TPO/PTR),
  `ptm_methyl.rtp` (MLZ/MLY/M3L), `gace.rtp` (GACE), `tym.rtp` (TYM), `pmx_bonded.itp` (TYM GAFF).
  Atom types clean. GACE/TYM GAFF types appended at END of atomtypes.atp/ffnonbonded.itp (marked).
- **`--reverse` feature**: committed (139922e0) — `pmx mutate --reverse` builds product-as-state-A.

## Key facts learned (don't re-derive)
- **Engine PTM codes** (`covalent/ptm.py`): SER→SP1/SP2, THR→TP1/TP2, TYR→YP1/YP2; LYS→me1/me2/me3.
  **pTyr is YP1/YP2** (I wasted time searching PDB name PTR — it's not the engine target).
- **Amber charges are shared across all variants** (LYS/ASP identical in 14sb/99sb/19sb/bsc1) — so
  ONE amber PTM set ports to all amber variants; variants differ only in dihedrals/bonded.
- **ECC scaling** (`-q-ecc` vs `-q`): holds peptide unit {N,H,C,O} FIXED, scales CA+sidechain so
  residue **net = ±0.75** (neutral residues untouched). Per-residue factor (LYS 0.776, ASP 0.72).
  **Reverse = rescale CA+sidechain so net returns to the integer** (helper `/tmp/reverse_ecc.py`,
  should be moved into the clean dir). Only methyl-Lys needed reversing; phospho had non-ECC source;
  GACE is neutral (unscaled).
- **Non-ECC amber phospho source** = `-q-mut` variants (SP1..YP2 at integer −1/−2), canonical.
- Canonical amber SEP/TPO use **standard amber types** (no phosaa add-on); phospho_2025 uses a
  different phosaa (P=1.38 vs canonical 1.323) — we chose **canonical**.

## DONE — me1→MLZ convergence (2026-08-11)
All four call sites now agree on **MLZ/MLY/M3L**:
1. ✅ `covalent/ptm.py`: `"me1": "MLZ"` (was MLYS).
2. ✅ charmm clean `ptm_methyl.rtp`: already `[ MLZ ]` (+ amber clean RTPs). No MLYS left in engine/clean.
3. ✅ `mutff_clean/smoke_test.sh` CASES: `LYS:me1:MLZ`.
4. ✅ `tests/test_ptm.py::test_methylation_degrees`: expects `("MLZ", 0)`.
**Golden guard re-run**: `131 passed / 1 failed` — exactly the baseline. The lone failure is the
PRE-EXISTING `test_golden.py::test_golden_hybrids_unchanged` term-del drift (cterm_del + nterm_del
emit no `hyb.top`); unrelated to me1. The me1 edit added ZERO new golden drift (methyl cases aren't
in the golden MATRIX, so all point/proton/nsaa cases stayed `ok`).
**Harness invocation** (record for next time): put devpmx's pmx on PATH but force PYTHONPATH so the
`pmx gentop` subprocess resolves to Documents/pmx/src (NOT off-limits GitHub/pmx):
```
export PATH="/Users/carter/miniforge3/envs/devpmx/bin:$PATH"
export PYTHONPATH=/Users/carter/Documents/pmx/src
export GMXLIB=/Users/carter/Documents/pmx/src/pmx/data/mutff
python3 tests/golden/generate_golden.py --check   # then: python3 -m pytest -q
```
Stale `MLYS` strings remain only in `_smoke/` build artifacts + prose in this file / CLEAN_FF_PLAN.md.

## Smoke test — `mutff_clean/smoke_test.sh <ff>`
Builds capped peptide → `pmx mutate` to PTM → pdb2gmx → gentop → editconf → grompp; reports missing
params. Matrix = {SER,THR,TYR}×{p1,p2} + LYS×{me1,me2,me3} on all three clean FFs.

### amber99sb-star-ildnp-q-mut added (2026-09-24) — 11/11
The 'q' variant of amber99 (charge-modified backbone, Best/Hummer + ILDN**P**
proline fix — the `p` in ildnp).  Built like the other clean -q FFs but from the
**amber99 family only** (amber 19 vs 99 differ: amber99 is NOT CMAP-based, and its
phospho is the OLDER paper — YP1 CZ=CA/OH=OS/P=+1.093, standard amber types, NOT
phosaa19SB).  Recipe: copy OG `amber99sb-star-ildnp-q-mut.ff` base (phospho + TYM
already baked in, GAFF bonded inline) → strip only S2SE.rtp/S2SP.rtp/genlib.log
(pre-built hybrids + log; keep gbsa.itp & aminoacids.vsd — forcefield.itp #includes
gbsa) → add amber99-family `ptm_methyl.rtp` + `gace.rtp` (CG=CT) from the amber99-std
clean.  Smoke 11/11; charges correct (SP1 dQ=-1, MLZ dQ=0, TYM dQ=-1, integer
endpoints).

### charmm36m-ecc-mut added (2026-09-24) — 9/9  [ECC re-layer #1]
First ECC clean FF.  **Verified CHARMM ECC rule: uniform sidechain × 0.75** (every
charged sidechain atom ×0.750, backbone untouched — measured std-LYS vs OG-ecc-LYS).
Recipe: copy OG `charmm36m-ecc-mut` base (phospho already baked in *and* ecc-scaled)
→ strip cruft.  Then the methyl fix (OG ecc's methyl was never given the clean-std
treatment — it had the OLD name `MLYS` + OLD protein typing CT2/NH3/CT3 that needs a
nonexistent angle): removed MLYS/MLY/M3L from the ecc base merged.rtp, and added an
**ecc-scaled `ptm_methyl.rtp` overlay** = clean-std CGenFF-typed methyl with sidechain
charges ×0.75 (MLZ/MLY/M3L net +0.75).  Also ported the 2 di-methyl Nz junction terms
missing from the ecc ffbonded (`CG334-NG3P1-CG334` angle + `HGA3-CG334-NG3P1-CG334`
dihedral).  Charges ecc-consistent: SP1 dQ=-0.75, SP2 dQ=-1.5, MLZ/MLY dQ=0.
**Charged termini already correct** (ecc tdb: N-term NH3+ nets +0.75, C-term COO- -0.75
→ N-term Lys = +1.5), inherited free from the OG ecc base's tdb.
TODO: amber -q-ecc (×3) next, same recipe (ecc-scaled methyl overlay + neutral gace).

### TYM quinoid Cz-Oeta = 0.130 corrected across ALL 8 FFs (2026-09-25)
The clean set shipped an UNCORRECTED tyrosinate everywhere (grompp-clean but wrong):
charmm-ecc had stock `CG2R61` (not the quinoid), charmm-std had no TYM, all 6 amber had the
stock `c-o` 0.1219 bond. Authoritative source = `~/projects/2026_amber_trysoine/ff_qmTYM`
(QM ωB97X-3c p-cresolate, Cz-Oeta=0.130). pmx hybrid bonded lookup is TYPE-based
(`get_bond_param(a.type)`), so the fix must live in the atom TYPE, not a global bond override.
- **Amber (6):** new dedicated quinoid Cz type **`cq`** (LJ clone of `c`); ONLY `cq-o`=0.130,
  the other 7 terms (1 bond cc-cq, 4 angles, 1 dihedral wildcard X-cc-cq-X, 1 improper cc-cc-cq-o)
  are verbatim `c`-clones. Global `c-o` stays 0.1219 → **no xFP caveat** (unlike ff_qmTYM's global
  override). Per-FF kb respected (amber19sbmut cq-o kb=533627.4). Patcher: scratchpad/patch_amber_cq.py.
- **CHARMM (2):** `CG2R6Q` via `ff_qmTYM/patch_charmm_CG2R6Q.py` (charmm-ecc: flipped stale
  merged.rtp TYM CZ CG2R61->CG2R6Q -- the patch's own merged.rtp edit no-ops on indented sub-headers,
  did it directly; charmm-std: backfilled TYM as a `tym.rtp` overlay from ff_qmTYM + ported 18 CGenFF
  aromatic-junction bonded terms clean-std's reduced base lacked).
- Hybrid Cz-Oeta B-state = **0.130** in all (amber `0.1364->0.130`, charmm `0.1411->0.130`) --
  byte-identical to ff_qmTYM's README verification. Charges intact (TYM net -1.0 std / -0.75 ecc).
- Smoke gate now runs TYM for charmm too (`amber*|charmm*`). All 8 FFs quinoid: amber CZ=cq, charmm CZ=CG2R6Q.

### FULL MATRIX RESULT: **86 / 86 grompp-clean** ✅  (6 amber×11 + charmm×10 + charmm-ecc×10; TYM now tested on charmm)
- Base 9 PTMs (phospho ×6 + methyl ×3) on all 6 FFs; **+ TYM + GACE on the 5 amber FFs**.
- amber99sb-star-ildn-mut 11/11 · amber14sbmut 11/11 · **amber14sb-q-mut 11/11** · amber19sbmut 11/11 ·
  amber19sb-q-mut 11/11 · charmm36m-mut 9/9 (TYM/GACE are amber GAFF residues; n/a to charmm).
- **amber14sb-q-mut** (the -q charge variant, 2026-08-16) built like amber19sb-q-mut: OG -q base (phospho
  phosaa14SB + TYM already baked in) minus genlib cruft; + ptm_methyl.rtp + gace.rtp (CG=2C).

### ACE/NME cap naming STANDARDIZED (2026-08-16)
The newer amber builds (amber14sb-q-mut, amber19sbmut, amber19sb-q-mut) named the ACE/NME cap methyl
H's `H1/H2/H3` (rtp + hdb base `H`), while the older ffs, the template, and pmx's `rename_atoms` all
use `HH31/HH32/HH33`.  That mismatch forced the `aminoacids.arn` cap map + an engine alias.  **Fix:
renamed those 3 ffs to `HH31/HH32/HH33`** in `aminoacids.rtp` (atoms+bonds) and `aminoacids.hdb` (base
`H`→`HH3`); the obsolete ACE/NME `.arn` cap entries were removed.  ⚠️ GOTCHA: NME's methyl H *type* is
the amber `H1` atom type (ACE's is `HC`) -- a naive `H1→HH31` rename clobbers that type column; rename
the atom-NAME field only.
- Kept as defensive robustness (now dormant for these ffs, but handles user PDBs of either
  convention): `alchemy.mutate` cap-routing + `nsaa._atom_name_aliases` (H1↔HH31).
- **smoke_test.sh terminus input** is auto-detected from `aminoacids.n.tdb`: a residue-specific tdb
  (has a `[ None ]` block -> amber*-q-mut) gives ACE a short `-ter` menu (feed `1\n2\n`); an empty tdb
  (amber14sbmut/amber99/amber19sbmut/charmm) uses GROMACS's generic menu (`3\n4\n`).  This is the real
  signal -- it's NOT tied to the ACE naming (amber19sbmut has HH31-named ACE but an empty tdb -> 3\n4\n).

### GACE is NOT a modified ACE (corrected 2026-08-16)
GACE is the GSH **γ-glutamyl residue** (a whole glutamate: free α-ammonium N1 + free α-carboxylate C1,
linking to the next residue through its **γ-carbonyl** off CG) — chemically unrelated to the acetyl
cap ACE.  An earlier pass wrongly built it by mutating `ACE→GACE` (grafting a glutamate onto ACE's 3
atoms — garbage structure, only grompp-clean by luck).  **Now tested as a real residue**: build the
GSH backbone **GACE-CYS-GLY-NME** from the RTP (`_mk_gace_pep.py`) and run pdb2gmx→grompp (no
mutation/gentop).  C-terminus is NME-capped so the build doesn't depend on each ff's C-terminal
carboxylate naming; the GACE γ-linkage to CYS is still exercised.  Consequently:
- **Reverted** the two engine workarounds added only for the bogus ACE→GACE: the `alchemy.mutate`
  cap-routing and `nsaa._atom_name_aliases` (golden still 131/1 after revert).
- **Real bug fixed**: the topic RTPs (gace/tym/ptm_methyl/ptm_phospho) lacked a `[ bondedtypes ]`
  directive, so their residues were unusable *directly* in a molecule with base residues ("different
  'all dihedrals' type" error) — the hybrid path had masked it.  Prepended each amber ff's base
  `[ bondedtypes ]` to its topic RTPs.
- TYM stays a genuine mutation case (`TYR→TYM`, a real tyrosinate modification).

### TYM + GACE added (2026-08-16) — "use amber14sb as a guide"
TYM (GAFF tyrosinate, net −1) and GACE (GSH γ-glutamyl N-cap, net-neutral) had their residues in the
clean FFs but their **GAFF bonded was un-wired** (`pmx_bonded.itp` was never `#include`d and was
incomplete). Fixed by merging the GAFF bonded **inline into ffbonded.itp** (the way OG
`amber14sb-q-mut` / `amber14sb-q-ecc-mut` ship it), per-family:
- **GACE typing is per-family**: amber99 uses `CG=CT` (+`CT-c3` junction); amber14/19 use `CG=2C`
  (+`2C-c3`). Our clean copies were all amber99-style `CG=CT` → replaced amber14/19 GACE with the
  `CG=2C` residue from the matching `*-q-ecc-mut` guide (GACE is net-neutral so ecc = non-ecc). TYM is
  net −1 so it MUST come non-ecc (the clean TYMs already were, sum −1.000).
- GAFF bonded merged via `/tmp/merge_gaff.py` (lowercase-type filter): amber99 from
  `amber99sb-star-ildnp-q-ecc-mut`, amber14 from `amber14sb-q-mut`+`-q-ecc-mut`, amber19 from
  `amber19sb-q-ecc-mut`. GAFF atom types (c3/cd/cc/nz/…) already present.
- **Smoke harness**: cases now carry a resid (`buildWT:code:target[:resid]`); `pmx mutate` gains
  `--generative`; amber FFs get `"TYR:TYM:TYM:3"` + `"ALA:GACE:GACE:1"` (GACE is an N-cap → resid 1).

### Two ENGINE fixes for GACE (both golden-safe, 131/1 held)
1. **`alchemy.mutate`** — a non-standard residue (a cap like ACE, `moltype=='unknown'`) now routes to
   the generative engine *before* `get_mtp_file` (which only knows protein/dna/rna and raised). Lets
   `pmx mutate --generative --script '1 GACE'` build ACE→GACE via the CLI.
2. **`nsaa._assign_rtp_entries`** — atom-name alias fallback (`_atom_name_aliases`): a cap methyl H is
   `H1/H2/H3` in the ff19SB rtp but `HH31/HH32/HH33` in the model (pmx/PDB convention); without this
   the engine left ACE's H's untyped in the ACE→GACE hybrid on amber19. (amber19 also keeps its
   `aminoacids.arn` HH3x→Hx map for pdb2gmx on the *non*-GACE cases, where ACE stays HH31.)

### amber19 build (2026-08-15) — ff19SB is CMAP-based; two new clean FFs
Built **both** `amber19sbmut` (standard ff19SB charges) and `amber19sb-q-mut` (the `-q` variant — a
distinct charge model for LYS/ASP/GLU/ARG; neutral residues identical). Decision: **data-only first,
CMAP-on-hybrids deferred** (the generative engine drops CMAP for non-charmm — `mutate.py:237`
`if bCharmm else []`; generated amber19 hybrids therefore carry no CMAP → still grompp-clean, but the
mutated residue lacks its ff19SB backbone correction; the shipped phospho residues themselves also
carry no `[ cmap ]`, so nothing dangles).
- **Phospho** = **phosaa19SB** = phosaa14SB with α-carbon type **CX→XC** (SI: `frcmod.phosaa19SB` +
  `phosaa19SB.lib` in `~/Downloads/phosaaSI/`). Source = `amber19sb-q-mut.ff` (GROMACS port).
  - `amber19sb-q-mut` clean: copied OG base (already has phospho + phosaa19SB types/bonded + TYM);
    **stripped** the non-generative cruft that collides with the engine (`S2*.rtp` pre-built hybrids,
    `add_mutres_cmap.py`, `cmap.ioldtp`, `genlib.log`); added `ptm_methyl.rtp` + `gace.rtp`.
  - `amber19sbmut` clean: OG standard base + **imported** phosaa19SB (6 residues + 7 types
    CG/OQ/OR/OT/OV/OX/OZ + nonbonded + 19 bond/50 angle/124 dih) via `/tmp/port_phosaa19.py`; added
    methyl + gace; **also imported GAFF atom types** (c3/cf/hc/hn/hx/nc/nd/ne/nf/nz) the standard base
    lacked but `gace.rtp` needs (pdb2gmx validates all rtp types at load → global fail without them).
- **Two harness/compat fixes (needed for any ff19SB):**
  1. **`aminoacids.arn`** — amber19 caps name methyl H's `H1/H2/H3` but the template PDB (and ff14SB)
     uses `HH31/HH32/HH33`; amber19 shipped no `.arn`. Added one (copied ff14SB's N-term H→H1 map +
     `ACE/NME HH3x→Hx`). Without it pdb2gmx: "Atom HH31 in ACE not found".
  2. **`smoke_test.sh` terminus input** — ff19SB's ACE/NME `-ter` menus are shorter, so the hardcoded
     `3\n4\n` is out-of-range → pdb2gmx re-prompts and **hangs** on the closed pipe. Now per-ff:
     `case "$FF" in amber19*) TER='1\n2\n'` (the "None/None" cap indices) `;; *) '3\n4\n'`. (Caps MUST
     be None — can't apply NH3+ to an acetyl cap.)
- amber19 pdb2gmx/grompp are noticeably slower (large CMAP) — one ff ≈ several minutes.

### What got it there (fixes applied this session)
- `_mk_wt_peptide.py` — NEW. Builds clean single-state WT capped peptide by grafting a
  `TargetResidue.from_rtp(WT)` residue onto the CYM template backbone (backbone-only Kabsch fit).
  Closes the harness gap. Signature: `_mk_wt_peptide.py <tmpl> <WT> <out.pdb> <ff>`.
- **smoke_test.sh** — (a) seeds the generated hybrid name into `$GMXLIB/residuetypes.dat` after
  mutate (the engine registers hybrids in the *packaged* mutff residuetypes.dat, not the GMXLIB in
  use — that's why THR/TYR read as "type Other"); (b) strips `nsaa.{rtp,mtp,arn}` per run for
  reproducibility (else stale/colliding hybrids leak via pdb2gmx's `*.rtp` glob).
- **ENGINE `covalent/ptm.py`**: `_WT_ALIAS = {"LYP":"LYS"}` — the charged lysine is renamed LYP
  before PTM resolve; it IS the methylation substrate, so resolve it. (test added)
- **ENGINE `nsaa.py`**: `_RTP_PROTON_ALIAS = {"LYP":"LYS"}` in `_assign_rtp_entries` — amber's
  charged-lysine RTP name is LYS, so a `me*` on LYP finds params. (Real amber bug, not smoke-only.)
- **amber PHOSPHO PROVENANCE (important, corrected 2026-08-11)** — the two amber clean FFs must use
  DIFFERENT phosphate papers, per Carter:
  - **amber99sb-star-ildn-mut** → older paper (`amber99sb-star-ildnp-q-mut.ff`). Correct as-is:
    CZ=CA, OH=OS, P=+1.093, pTyr aromatic bond CA-OS. I appended that pTyr CA-OS block here.
  - **amber14sbmut** → **phosaa14SB** (Homeyer/Steinbrecher, newer paper; SI = `frcmod.phosaa14SB` +
    `phosaa14SB.lib` in `~/Downloads/phosaaSI/`). The clean ff had WRONGLY inherited the amber99
    params (byte-identical) — NOT evident in the files. **Now rebuilt** from `amber14sb-q-mut.ff`
    (verified = GROMACS port of frcmod.phosaa14SB: `C-OR`=450/1.323→0.1323/376560, `CA-C-OR-P`
    n1=17.975×4.184=75.209 both match). Migration (`/tmp/port_phosaa.py`): replaced ptm_phospho.rtp
    with phosaa14SB SP1–YP2 (types C/OR/OQ/OP mono, OT/OV/OX/OZ di, CG Tyr-Cγ; P=+1.393); added 7
    atom types + nonbonded (FEP-tuned LJ) + phosaa bonded (19 bond / 50 angle / 124 dihedral); removed
    the vestigial amber99 CA-OS block. amber14sbmut re-validated **9/9**. Net charges integer (−1/−2),
    dianions use OT/OV (pTyr²⁻) & OX/OZ (pSer/pThr²⁻).
  - **amber19** (not in clean tree yet): use **phosaa19SB** (`frcmod.phosaa19SB` attached; ≈14SB) via
    `amber19sb-q-mut.ff` the same way when its clean ff is built.
- **charmm dedup**: removed pdb2gmx `*.rtp`-glob collisions — `MLZ` dup (ptm_extra vs ptm_methyl),
  `NME` (ptm_misc vs base merged.rtp), and `TP1/TP2/TPO` (ptm_phospho vs base merged.rtp). The OG
  charmm `merged.rtp` NATIVELY ships TP1/TP2/TPO+NME; only SEP/PTR/YP1/YP2 are pmx additions.
- **charmm `ptm_methyl.rtp` MLZ**: re-typed to CGenFF (NG3P2/CG324/CG334, methyl carbon CM) to match
  MLY/M3L. The old standard-protein-typed MLZ (CT2/NH3/CT3) needed a non-existent CT2-NH3-CT3 angle.
  NOTE: our atomtyping is CORRECT — it matches canonical CHARMM (`aminoacids.rtp` ships MLZ/MLY/M3L
  with exactly CB=CG321, CE=CG324, NZ=NG3P2/P1/P0, methyl=CG334).
- **charmm ffbonded.itp — methyl-Lys junction block appended** (the fix that closed 24→27). The base
  `charmm36m-mut.ff` in the tree is a *reduced* CHARMM: MLZ/MLY/M3L are in `merged.rtp`, but the
  CGenFF-modified-residue junction bonded was never pulled in. Ported verbatim from the full CHARMM
  distribution on the crucialP310 volume (`charmm36-feb2026_cgenff-5.0.ff/ffbonded.itp`). Exactly the
  terms the hybrids need (enumerated by instrumenting `BondedParser.get_*_param`, NOT by grep which
  over-pulled 40 KB): backbone↔sidechain seam **bond** CT1-CG321; **angles** CG321-CT1-{C,HB1,NH1},
  CT1-CG321-{CG321,HGA2}; di-methyl Nz **angle** CG334-NG3P1-CG334; di-methyl Nz **dihedral**
  HGA3-CG334-NG3P1-CG334. (mono MLZ + tri M3L needed only the seam terms; di MLY needed the two Nz
  terms too.) These are canonical CHARMM values (provenance `ALY`/`nc3 fylin` etc.), not invented.

### Latent bug (found, not yet fixed) — hybrid-name collision
MLZ & MLY both truncate to hybrid name `K2ML` (`(K+2+MLZ)[:4]`), so they overwrite each other in
`nsaa.*`. Harmless in the sequential smoke (each case rebuilds K2ML before its own gentop, and the
harness now strips nsaa.* per run) but WILL corrupt a system that contains both a mono- and a
di-methyl-Lys. Fix in the engine (`build_hybrid_residue` naming) — see TODO 2.

## TODO (ordered)
1. **Hybrid-name collision** MLZ/MLY→K2ML: give methyl-Lys targets distinct 4-char hybrid names so
   both can coexist in one system.
2. TYM / GACE cases are NOT yet in the smoke matrix — add them (bonded coverage for GAFF types).
3. **ECC re-layer**: build `charmm36m-ecc-mut` clean = same overlay, ECC-scale charged side chains.
   Amber ECC likewise. (Note: amber clean methyl-Lys RTP net charge reads 0.75 — ECC-scaled — verify
   whether the non-ECC clean amber should be integer +1.)
4. Move `/tmp/reverse_ecc.py` + `/tmp/extract_rtp.py` (comment-tolerant) into the clean dir.

## Guardrails
- OG `mutff/` and `GitHub/pmx` are read-only references. Work only in `mutff_clean/`.
- `PYTHONPATH=/Users/carter/Documents/pmx/src` (else off-limits GitHub/pmx).
- Golden baseline = 131 pass / 1 fail (term_del, pre-existing).

### amber -q-ecc trio built (2026-09-25) — 11/11 each ✓  [ECC re-layer COMPLETE]
amber14sb-q-ecc-mut, amber19sb-q-ecc-mut, amber99sb-star-ildnp-q-ecc-mut.
Recipe was trivial: OG -q-ecc bases already bake in phospho+TYM+**methyl+GACE** (all ecc-scaled,
LYS +0.75), so just: copy base -> strip S2*.rtp/genlib.log/nsaa.*/vsd (keep gbsa) -> run
scratchpad/patch_amber_cq.py (cq TYM fix, cq-o=0.130). No methyl/gace overlays needed (baked in).
All 11/11, SP1 dQ=-0.75, TYM bond 0.130. GOTCHA fixed: amber14/19 -q-ecc bases had OLD ACE/NME cap
naming (H1/H2/H3) never standardized like the -q siblings -> pdb2gmx "HH31 not found in ACE";
ported ACE+NME rtp+hdb from the -q sibling (caps are unscaled -> identical bar names). ACE cap naming
now uniform HH31 across all 11 FFs; no stale .arn. **Clean set = 11 FFs, matrix 119/119.**

### harness relocated (2026-09-25)
smoke_test.sh + _mk_wt_peptide.py + _mk_gace_pep.py MOVED to `tests/mutff_clean_smoke/`
(tooling, not FF content). It sets GMXLIB to this data dir itself. Run:
`bash tests/mutff_clean_smoke/smoke_test.sh <ff>`. See that dir's README.md.

### upstream alignment: +3 pristine base FFs + aux files (2026-09-27)
Clean set was missing 3 of upstream deGrootLab's 6 FFs. Pulled them in VERBATIM from
`upstream/develop` (git archive; 0-file-diff verified pristine): amber99sb-star-ildn-bsc1-mut,
amber99sb-star-ildn-dna-mut, charmm22star-mut. These are base FFs (nucleic + old CHARMM) with
mutres for STANDARD mutations only -- no PTM overlays, so the PTM smoke matrix doesn't apply
(they're deGroot-validated as-is). Also added the shared aux files (specbond/atommass/elements/
vdwradii/xlateat/spc216.gro) so the set is self-contained. amber .arn ACE/NME cap maps standardized.
**Clean set = 14 FFs = upstream-6 (pristine) + fork-8 (-q/-ecc/amber19/PTM, TYM-quinoid). Only the
fork-8 are PTM-capable (119/119 smoke).** Dropped only fork-only amber99sbdisp + des-amber.
