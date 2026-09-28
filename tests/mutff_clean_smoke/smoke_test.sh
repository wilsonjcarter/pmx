#!/usr/bin/env bash
# smoke_test.sh <ff-name>  — for each PTM: build ACE-ALA-<WT>-ALA-NME, pmx-mutate to the
# PTM target, pdb2gmx -> gentop -> grompp, and report missing FF params.
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
FF="${1:?usage: smoke_test.sh <ff-name e.g. amber99sb-star-ildn-mut>}"
export PYTHONPATH=/Users/carter/Documents/pmx/src
# This harness lives in tests/; the clean FFs live in the pmx data dir.
FFROOT="$(cd "$HERE/../../src/pmx/data/mutff_clean" && pwd)"
export GMXLIB="$FFROOT"
PY=/Users/carter/miniforge3/envs/devpmx/bin/python
pmx(){ "$PY" -m pmx.scripts.cli "$@"; }
TMPL=/Users/carter/Documents/pmx/examples/02_protonation_C2CM/input/ref_peptide.pdb  # ACE-ALA-CYM-ALA-NME
# Strip engine-accumulated hybrid artifacts so each matrix run is reproducible
# (the generative engine upserts hybrids into nsaa.{rtp,mtp,arn}; leftovers from
# a prior run -- or a name-colliding sibling -- otherwise leak in via pdb2gmx's
# *.rtp glob).  Same rationale as generate_golden.py's _strip_runtime_hybrids.
rm -f "$FFROOT/$FF.ff"/nsaa.rtp "$FFROOT/$FF.ff"/nsaa.mtp "$FFROOT/$FF.ff"/nsaa.arn
W="$HERE/_smoke/$FF"; rm -rf "$W"; mkdir -p "$W"; cd "$W"
# "buildWT:code:target[:resid]"  -- build ACE-ALA-<buildWT>-ALA-NME, then mutate
# residue <resid> (default 3) to <code> (a PTM alias p1/me1... OR a direct NSAA
# target name like TYM/GACE).  GACE is an N-terminal cap, so it mutates resid 1.
declare -a CASES=(
  "SER:p1:SP1" "SER:p2:SP2" "THR:p1:TP1" "THR:p2:TP2" "TYR:p1:YP1" "TYR:p2:YP2"
  "LYS:me1:MLZ" "LYS:me2:MLY" "LYS:me3:M3L"
)
# TYM (GAFF tyrosinate) is a real TYR->TYM modification -> mutation case.
# GACE is NOT a mutation target (it's the GSH gamma-glutamyl residue) -> tested
# separately below by building real glutathione (GACE-CYS-GLY).
case "$FF" in amber*|charmm*) CASES+=( "TYR:TYM:TYM:3" );; esac
pass=0; fail=0
for c in "${CASES[@]}"; do
  IFS=: read wt code tgt rid <<<"$c"; rid=${rid:-3}
  d="$W/${wt}_${tgt}"; mkdir -p "$d"; cd "$d"
  # build ACE-ALA-<wt>-ALA-NME : mutate the CYM template (resid 3) -> wt
  "$PY" "$HERE/_mk_wt_peptide.py" "$TMPL" "$wt" wt.pdb "$FF" >mkwt.log 2>&1   # clean WT peptide (no hybrid)
  # mutate residue <rid> -> target (generative; --generative also routes direct
  # NSAA targets and cap mutations like ACE->GACE that lack a pre-built hybrid)
  printf '%s %s\n' "$rid" "$code" > mut.txt
  pmx mutate --generative -f wt.pdb -o hyb.pdb -ff "$FF" --script mut.txt >mut.log 2>&1
  # The engine registers the generated hybrid in the *packaged* mutff
  # residuetypes.dat (hardcoded), not the GMXLIB in use.  When GMXLIB is this
  # clean tree, seed the name here so pdb2gmx sees the hybrid as a Protein
  # residue rather than "type Other".
  hyb=$(sed -n 's/.*hybrid \([A-Za-z0-9]*\) .*/\1/p' mut.log | head -1)
  if [ -n "$hyb" ] && ! grep -qw "$hyb" "$GMXLIB/residuetypes.dat"; then
    [ -n "$(tail -c1 "$GMXLIB/residuetypes.dat")" ] && echo >> "$GMXLIB/residuetypes.dat"
    printf '%-6s Protein\n' "$hyb" >> "$GMXLIB/residuetypes.dat"
  fi
  # ACE/NME caps must take the "None" terminus, but its menu index differs by ff
  # -- an out-of-range index makes pdb2gmx re-prompt and hang on the closed pipe.
  # ffs with a *residue-specific* aminoacids.n.tdb (a '[ None ]' block) give ACE
  # a short menu (None at 1/2); ffs with an empty n.tdb fall back to GROMACS's
  # generic menu (None at 3/4).  Key off the tdb -- ff-name-independent.
  if grep -q '^\[ None \]' "$HERE/$FF.ff/aminoacids.n.tdb" 2>/dev/null; then
    TER='1\n2\n'; else TER='3\n4\n'; fi
  printf '%b' "$TER" | gmx pdb2gmx -f hyb.pdb -o hg.pdb -p htop.top -ff "$FF" -water none -ter >p2.log 2>&1
  pmx gentop -p htop.top -o pmx.top -ff "$FF" >gt.log 2>&1
  printf 'integrator=steep\nnsteps=1\ncutoff-scheme=Verlet\nfree-energy=yes\ninit-lambda-state=0\nfep-lambdas=0 1\nsc-alpha=0.3\nsc-power=1\nsc-sigma=0.25\n' > em.mdp
  gmx editconf -f hg.pdb -o box.pdb -bt cubic -d 1.5 >ec.log 2>&1
  gmx grompp -f em.mdp -c box.pdb -p pmx.top -o em.tpr -maxwarn 10 >grompp.log 2>&1
  if [ -f em.tpr ]; then echo "  OK    $wt -> $tgt ($code)"; pass=$((pass+1))
  else
    err=$(grep -iE 'not found|no default|Fatal' */*.log grompp.log mut.log gt.log 2>/dev/null | grep -ivE 'warning' | head -1)
    echo "  FAIL  $wt -> $tgt ($code)  :: ${err:-see logs in $d}"; fail=$((fail+1))
  fi
done

# --- GACE tested as a REAL residue: build glutathione (GACE-CYS-GLY) ---------
# GACE is the GSH N-terminal gamma-glutamyl residue (linked to CYS via its
# gamma-carbonyl); it is NOT a modified ACE, so it is validated by building the
# GSH backbone GACE-CYS-GLY(-NME) from the RTP and running pdb2gmx -> grompp (no
# mutation/gentop).  The C-terminus is NME-capped so the build doesn't depend on
# each ff's C-terminal carboxylate naming; GACE's gamma-linkage is still exercised.
case "$FF" in amber*)
  d="$W/GSH_GACE"; mkdir -p "$d"; cd "$d"
  "$PY" "$HERE/_mk_gace_pep.py" "$FF" gsh.pdb >mk.log 2>&1
  # GACE/NME cap terminus menu indices vary by ff -- brute-force the pair,
  # guarding each try against a hang from an out-of-range index.
  for T in '3\n2\n' '1\n2\n' '0\n2\n' '3\n0\n' '0\n0\n'; do
    rm -f g.top
    ( printf '%b' "$T" | gmx pdb2gmx -f gsh.pdb -o g.gro -p g.top -ff "$FF" -water none -ter >p2.log 2>&1 ) &
    pid=$!; ( sleep 30 && kill -9 $pid 2>/dev/null ) & wd=$!
    wait $pid 2>/dev/null; kill $wd 2>/dev/null; wait $wd 2>/dev/null
    [ -s g.top ] && break
  done
  if [ -s g.top ]; then
    gmx editconf -f g.gro -o gb.gro -bt cubic -d 1.5 >ec.log 2>&1
    printf 'integrator=steep\nnsteps=1\ncutoff-scheme=Verlet\n' > gem.mdp
    gmx grompp -f gem.mdp -c gb.gro -p g.top -o g.tpr -maxwarn 20 >grompp.log 2>&1
  fi
  if [ -f g.tpr ]; then echo "  OK    GACE  (GSH backbone GACE-CYS-GLY-NME)"; pass=$((pass+1))
  else echo "  FAIL  GACE  (GSH)  :: $(grep -hiE 'fatal|not found|no default' "$d"/p2.log "$d"/grompp.log 2>/dev/null | grep -iv warning | head -1)"; fail=$((fail+1)); fi
  ;;
esac
echo "== $FF: $pass ok / $fail fail =="
