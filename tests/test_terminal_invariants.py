"""Invariant + contract tests for terminal edits (deletion, amidation, terminal mutation).

Complements ``tests/golden``: golden pins the *exact bytes* of a few CHARMM cases,
this pins *properties* that must hold for any force field, so an AMBER regression or a
new residue class fails here even with no golden fixture.

Three groups:

* force-field-agnostic contract -- ``_rtp_terminal`` must read the FF's own ``.tdb``
  and raise, never silently fall back to another force field's parameters.
* ``.tdb`` parsing -- both five-token ``[ replace ]`` forms.
* generated-topology invariants -- run only when ``gmx`` is available, each in an
  isolated force-field copy so the tracked FF is never written to.
"""
import collections
import os
import re
import shutil
import subprocess
import sys
import tempfile

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))

from pmx.nsaa import (_rtp_terminal, _apply_tdb_variant, _tdb_path, _parse_tdb,
                      _tdb_variant_name)
from pmx.gmx import set_gmxlib

PEPTIDE = os.path.join(HERE, "golden", "_peptide.pdb")
needs_gmx = pytest.mark.skipif(shutil.which("gmx") is None, reason="gmx not on PATH")


@pytest.fixture(scope="module", autouse=True)
def _gmxlib():
    set_gmxlib()
    return os.environ["GMXLIB"]


# ---------------------------------------------------------------- contract ----

def test_rtp_terminal_requires_the_force_field():
    """No ff/resname must raise, not silently use a hardcoded patch."""
    with pytest.raises(ValueError, match="needs the force field"):
        _rtp_terminal({"atoms": [["N", "N", 0.0, 1]], "bonds": []}, "N")


def test_rtp_terminal_raises_when_neither_convention_applies():
    """A terminus is buildable two ways, and pmx tries both before failing:

    * a patch from the force field's ``.n/.c.tdb`` (charmm*, amber*-q-*), tried
      FIRST because the -q- ports ship both and the tdb is their validated path;
    * an explicit terminal RTP entry, ``NALA`` / ``CALA`` (the plain amber ports,
      whose 10-byte placeholder .tdb files are deliberate, not a gap).

    It raises only when NEITHER exists -- here a PTM residue, for which no port
    defines ``NMLZ``/``CMLZ`` and no tdb variant applies.
    """
    with pytest.raises(ValueError, match="defines no 'N' variant"):
        _rtp_terminal({"atoms": [["N", "N", 0.0, 1]], "bonds": []}, "N",
                      ff="amber14sbmut", resname="MLZ")


def test_placeholder_tdb_falls_back_to_the_explicit_entry():
    """A placeholder tdb is NOT a gap: amber14sbmut carries NALA/CALA instead."""
    entr = _rtp_terminal({"atoms": [["N", "N", 0.0, 1]], "bonds": []}, "N",
                         ff="amber14sbmut", resname="ALA")
    names = [a[0] for a in entr["atoms"]]
    assert "N" in names and len(names) > 5, names      # the real NALA entry
    assert any(n in names for n in ("H1", "H2", "H3")), names


@pytest.mark.parametrize("ff", ["charmm36m-mut", "amber14sb-q-mut"])
@pytest.mark.parametrize("term,expect", [("N", ("NH3+", "NH2+")), ("C", ("COO-",)),
                                         ("C_amide", ("CT2",))])
def test_conventional_variant_names_are_present(ff, term, expect):
    """Both supported families must use the conventional GROMACS variant names --
    that convention is the agnostic contract (we do not scan tdb semantics)."""
    side = "N" if term == "N" else "C"
    variants = _parse_tdb(_tdb_path(ff, side))
    name = _tdb_variant_name(variants, term, "ALA")
    assert name is not None, f"{ff} has no {term} variant; present: {sorted(variants)}"
    assert any(e in name for e in expect), (name, expect)


# ------------------------------------------------------------- tdb parsing ----

def test_tdb_replace_accepts_both_five_token_forms():
    """GROMACS has two five-token [ replace ] rows, told apart only by whether the
    third token is numeric.  Reading the no-rename form as a rename put a MASS
    string in the atomtype slot ('12.011'), which downstream became element '1'
    and crashed the RDKit embed.
        old NEW  type mass q [cgnr]   -- amber COO-: O OC1 O2     16.00  -0.7908
        old type mass q     [cgnr]    -- amber CT2 : C C   12.011  0.5973  1
    """
    entr = {"atoms": [["C", "XX", 0.0, 1], ["O", "YY", 0.0, 2]],
            "bonds": [], "diheds": [], "improps": [], "cmap": []}
    rename = [("replace", [["O", "OC1", "O2", "16.00", "-0.7908"]])]
    got = _apply_tdb_variant(entr, rename, "C")
    assert ["OC1", "O2"] == [got["atoms"][1][0], got["atoms"][1][1]]

    norename = [("replace", [["C", "C", "12.011", "0.5973", "1"]])]
    got = _apply_tdb_variant(entr, norename, "C")
    row = next(a for a in got["atoms"] if a[0] == "C")
    assert row[1] == "C", f"atomtype must be 'C', not the mass string: {row}"
    assert abs(float(row[2]) - 0.5973) < 1e-9, row


# --------------------------------------------------- generated topologies ----

_RTP_SECTIONS = {"bondedtypes", "atoms", "bonds", "angles", "dihedrals",
                 "impropers", "cmap", "exclusions"}


def _block_names(path):
    if not os.path.isfile(path):
        return set()
    txt = open(path, errors="replace").read()
    return set(re.findall(r"^\[\s*(\S+)\s*\]", txt, re.M)) - _RTP_SECTIONS


def _dedup_vs_library(ffdir):
    """Drop from the shipped mutres/mutres_term libraries any residue block the
    freshly generated nsaa.* now owns (mirrors terminal/_dedup_vs_mutres.py)."""
    for ext in ("rtp", "mtp"):
        owned = _block_names(os.path.join(ffdir, "nsaa." + ext))
        if not owned:
            continue
        for lib in ("mutres_term." + ext, "mutres." + ext):
            path = os.path.join(ffdir, lib)
            if not os.path.isfile(path):
                continue
            out, skip, removed = [], False, 0
            for l in open(path, errors="replace"):
                m = re.match(r"^\[\s*(\S+)\s*\]", l)
                if m:
                    skip = m.group(1) in owned
                    if skip:
                        removed += 1
                        continue
                if not skip:
                    out.append(l)
            if removed:
                open(path, "w").writelines(out)


def _run(cmd):
    """Run *cmd*, failing the test with its stderr rather than a bare exit code."""
    r = subprocess.run(cmd, stdin=subprocess.DEVNULL, capture_output=True, text=True)
    if r.returncode != 0:
        tail = "\n".join((r.stderr or r.stdout or "").strip().splitlines()[-12:])
        raise AssertionError("%s failed:\n%s" % (" ".join(cmd[:4]), tail))
    return r


def _isolated_ff(ff, tmpdir):
    """Copy one FF + the top-level GMXLIB files so a build cannot touch the tree."""
    src = os.environ["GMXLIB"]
    root = os.path.join(tmpdir, "ff")
    os.makedirs(root, exist_ok=True)
    for f in os.listdir(src):
        p = os.path.join(src, f)
        if os.path.isfile(p):
            shutil.copy(p, root)
    shutil.copytree(os.path.join(src, ff + ".ff"), os.path.join(root, ff + ".ff"))
    for stale in ("nsaa.rtp", "nsaa.mtp", "nsaa.arn"):
        p = os.path.join(root, ff + ".ff", stale)
        if os.path.isfile(p):
            os.remove(p)
    return root


def _truncated_peptide(path, first, last):
    """Write residues [first, last] of the golden peptide (1-based) to *path*.

    The golden peptide is ALA VAL ASP GLU LEU ILE LYS HSD CYS.  Residues 8-9 use
    CHARMM-specific names (``HSD``) that AMBER's RTP does not define, so the
    force-field-agnostic tests use the 1..7 window only.
    """
    with open(path, "w") as out:
        for l in open(PEPTIDE):
            if l[:4] == "ATOM" and first <= int(l[22:26]) <= last:
                out.write(l)
        out.write("END\n")


def _build(ff, op, tmpdir, first=1, last=7, target=None):
    """Apply *op* and run pdb2gmx+gentop.  Returns the topology text."""
    from pmx.model import Model
    root = _isolated_ff(ff, tmpdir)
    cwd, saved = os.getcwd(), os.environ.get("GMXLIB")
    try:
        os.chdir(tmpdir)
        os.environ["GMXLIB"] = root
        _truncated_peptide("start.pdb", first, last)
        _run(["gmx", "pdb2gmx", "-f", "start.pdb", "-o", "wt.pdb",
              "-p", "wt.top", "-ff", ff, "-water", "tip3p", "-ignh"])
        m = Model("wt.pdb", bPDBTER=True, rename_atoms=True)
        res = m.residues
        if op == "delN":
            from pmx.covalent import delete_terminal
            delete_terminal(m, res[0].id, ff, inplace=True)
        elif op == "delC":
            from pmx.covalent import delete_terminal
            delete_terminal(m, res[-1].id, ff, inplace=True)
        elif op == "amidate":
            from pmx.covalent import amidate_cterminus
            amidate_cterminus(m, res[-1].id, ff, inplace=True)
        elif op == "mutate":
            from pmx.alchemy import mutate
            mutate(m, res[-1].id, target, ff, inplace=True)
        else:
            raise ValueError(op)
        m.write("mut.pdb")
        # The generative engine writes e.g. VdeN into nsaa.rtp while CHARMM also
        # SHIPS VdeN in mutres_term.rtp -> pdb2gmx aborts with "found rtp entries
        # for X in both".  The production build scripts resolve this with
        # terminal/_dedup_vs_mutres.py: keep the fresh generative entry, drop the
        # stale library copy.  Do the same here.
        _dedup_vs_library(os.path.join(root, ff + ".ff"))
        _run(["gmx", "pdb2gmx", "-f", "mut.pdb", "-o", "gmx.pdb",
              "-p", "pmxtop.top", "-ff", ff, "-water", "tip3p"])
        _run(["pmx", "gentop", "-p", "pmxtop.top", "-o", "fe.top", "-ff", ff])
        return open("fe.top").read(), os.path.join(root, ff + ".ff")
    finally:
        os.chdir(cwd)
        if saved is not None:
            os.environ["GMXLIB"] = saved


def _atoms_and_bonded(top_text):
    cur, S = None, collections.defaultdict(list)
    for l in top_text.splitlines():
        s = l.split(";")[0].strip()
        if s.startswith("["):
            cur = s.strip("[] ").strip()
            continue
        if cur and s:
            S[cur].append(s.split())
    return S


def _ff_atomtypes(ffdir):
    T, cur = set(), None
    p = os.path.join(ffdir, "ffnonbonded.itp")
    for l in open(p, errors="replace"):
        s = l.split(";")[0].strip()
        if s.startswith("["):
            cur = s.strip("[] ").strip()
            continue
        if cur == "atomtypes" and s:
            T.add(s.split()[0])
    return T


@needs_gmx
@pytest.mark.parametrize("ff", ["charmm36m-mut", "amber14sb-q-mut"])
# The 1..7 peptide is ALA VAL ASP GLU LEU ILE LYS: its C-terminal residue is LYS
# (formal +1), so removing or replacing it CHANGES the net charge.  A terminal
# deletion is only charge-conserving when the deleted residue is neutral (which is
# why the degron ladder -- S/P/M/A -- needs no double-system-single-box).
@pytest.mark.parametrize("op,target,dq", [("delN", None, 0.0),    # deletes ALA1, neutral
                                          ("delC", None, -1.0),   # deletes LYS7, +1
                                          ("amidate", None, None),
                                          ("mutate", "ASP", -2.0),  # LYS +1 -> ASP -1
                                          ("mutate", "LYS", 0.0),   # LYS -> LYS
                                          ("mutate", "GLY", -1.0),  # LYS +1 -> GLY 0
                                          ("mutate", "PRO", -1.0)]) # LYS +1 -> PRO 0
def test_generated_topology_invariants(ff, op, target, dq, tmp_path):
    top, ffdir = _build(ff, op, str(tmp_path), target=target)
    S = _atoms_and_bonded(top)
    known = _ff_atomtypes(ffdir)
    rows = [r for r in S["atoms"] if r and r[0].isdigit()]
    assert rows, "no atoms parsed"

    qa = qb = 0.0
    for r in rows:
        tA, qA = r[1], float(r[6])
        hasB = len(r) >= 11
        tB, qB = (r[8], float(r[9])) if hasB else (tA, qA)
        qa += qA
        qb += qB
        for t, st in ((tA, "A"), (tB, "B")):
            # a blank or mass-like atomtype is the signature of an atom the RTP
            # entry never defined (the OC1/OC2 and OT1/OT2 bug class)
            assert t and not re.fullmatch(r"[-+0-9.eE]+", t), \
                f"{ff}/{op}: atom {r[4]} res{r[2]} state {st} has bogus atomtype {t!r}"
            assert t in known, f"{ff}/{op}: atomtype {t!r} absent from ffnonbonded"
        if hasB:
            # A mass change between states is LEGITIMATE: an element-changing map
            # is standard pmx practice, and upstream's own mutres library does it
            # 112 times for C<->O alone (S2C maps OG mass 16.000 -> SG 32.060).
            # What must never happen is the bug this guards -- an atom whose mass
            # was silently filled in from somewhere other than its own type -- so
            # require that a mass change comes WITH a type change.
            if abs(float(r[7]) - float(r[10])) >= 1e-3:
                assert tA != tB, (
                    f"{ff}/{op}: atom {r[4]} mass changed {r[7]} -> {r[10]} but "
                    f"its atomtype did not ({tA}); a silent mass change is the "
                    f"OC1/OC2-OT1/OT2 bug class, an element change is not")

    for q, st in ((qa, "A"), (qb, "B")):
        assert min(abs(q - round(q)), abs(q / 0.75 - round(q / 0.75))) < 2e-3, \
            f"{ff}/{op}: state {st} total charge {q:+.5f} neither integral nor 0.75-scaled"
    if dq is not None:
        assert abs((qb - qa) - dq) < 2e-3, f"{ff}/{op}: dq={qb - qa:+.4f}, expected {dq:+.1f}"

    # No bonded term may be listed twice with a NONZERO force constant.  Canonicalise
    # as the tuple or its REVERSAL -- never sorted(): i-j-k-l and j-i-k-l are
    # different dihedrals and sorting them together invents duplicates.
    for sec, nat in (("bonds", 2), ("angles", 3), ("dihedrals", 4)):
        seen = collections.Counter()
        for r in S.get(sec, []):
            if len(r) >= nat + 1 and r[0].isdigit():
                a = tuple(int(x) for x in r[:nat])
                seen[(min(a, a[::-1]), r[nat], tuple(r[nat + 1:]))] += 1
        for key, n in seen.items():
            if n < 2:
                continue
            nums = [x for x in key[2] if re.fullmatch(r"-?\d+(\.\d+)?", x)]
            ks = nums[1::3] if sec == "dihedrals" and len(nums) % 3 == 0 else nums[1::2]
            if ks and all(abs(float(k)) < 1e-12 for k in ks):
                continue          # duplicated k=0 term: 0+0=0, harmless
            pytest.fail(f"{ff}/{op}: {sec} {key[0]} listed {n}x with nonzero k "
                        f"(params {key[2]}) -- doubles the force constant")


@needs_gmx
@pytest.mark.parametrize("ff", ["charmm36m-mut", "amber14sb-q-mut"])
@pytest.mark.parametrize("op", ["delN", "delC"])
def test_deletion_leaving_one_residue_raises(ff, op, tmp_path):
    """A deletion that leaves a ONE-residue chain makes the survivor both N- and
    C-terminal, but it is typed from the internal RTP entry, so the terminus atoms
    it already carries (OT1/OT2 or H1/H2/H3) would stay untyped and surface as an
    unparseable three-field nsaa.rtp row inside pdb2gmx.  It must raise instead.
    """
    with pytest.raises(ValueError, match="no atomtype"):
        _build(ff, op, str(tmp_path), first=6, last=7)


@needs_gmx
@pytest.mark.parametrize("ff", ["charmm36m-mut", "amber14sb-q-mut"])
def test_deletion_state_b_equals_an_independent_build(ff, tmp_path):
    """State B of a terminal deletion must be exactly what pdb2gmx builds for the
    SHORTENED peptide: same (atomtype, charge) multiset outside the deleted residue.
    Atom NAMES differ by design (the hybrid keeps HV1/HV2 for appearing dummies).
    """
    hybrid, _ = _build(ff, "delN", str(tmp_path / "h"), first=1, last=7)
    # reference = plain pdb2gmx of residues 2..7, no alchemy
    ref = None
    root = _isolated_ff(ff, str(tmp_path / "ref"))
    cwd, saved = os.getcwd(), os.environ.get("GMXLIB")
    try:
        os.makedirs(str(tmp_path / "ref"), exist_ok=True)
        os.chdir(str(tmp_path / "ref"))
        os.environ["GMXLIB"] = root
        _truncated_peptide("start.pdb", 2, 7)
        _run(["gmx", "pdb2gmx", "-f", "start.pdb", "-o", "wt.pdb",
              "-p", "wt.top", "-ff", ff, "-water", "tip3p", "-ignh"])
        ref = open("wt.top").read()
    finally:
        os.chdir(cwd)
        if saved is not None:
            os.environ["GMXLIB"] = saved

    H, R = _atoms_and_bonded(hybrid), _atoms_and_bonded(ref)
    hb = collections.Counter()
    for r in H["atoms"]:
        if not (r and r[0].isdigit()) or r[2] == "1":      # skip the deleted residue
            continue
        tB = r[8] if len(r) >= 11 else r[1]
        qB = float(r[9]) if len(r) >= 11 else float(r[6])
        hb[(tB, round(qB, 4))] += 1
    wt = collections.Counter((r[1], round(float(r[6]), 4))
                             for r in R["atoms"] if r and r[0].isdigit())
    assert hb == wt, (f"{ff}: state B != independent build of the shortened peptide; "
                      f"only in B: {sorted((hb - wt).elements())[:6]}, "
                      f"only in WT: {sorted((wt - hb).elements())[:6]}")

@pytest.mark.parametrize("target", ["CYM", "SP1", "MLZ"])
def test_terminal_name_over_four_chars_is_refused(target):
    """A terminal hybrid name of 5 characters cannot round-trip through a PDB --
    it overflows resName (cols 18-21) into chainID and garbles resSeq:

        'ATOM  144  OT2 C2CMcX   9'  ->  resName='C2CM' chainID='c' resSeq='X'

    It happens whenever a side has a TWO-character library code (CYM=CM, TYM=YM,
    LYSN=KN, CYSO=CO, phospho P1/P2) or the target is a PTM/NSAA (``_target_code``
    is two characters).  pmx refuses rather than emitting a corrupt structure; see
    the memory note terminal-hybrid-naming-open-problem for the options.
    """
    from pmx.library import hybrid_name
    for term in ("N", "C"):
        with pytest.raises(ValueError, match="cannot round-trip through a PDB"):
            hybrid_name("CYS" if target == "CYM" else "SER", target, term=term)


def test_terminal_name_for_canonical_targets_fits_four_chars():
    """The canonical space must always fit, since that is the 99% case."""
    from pmx.library import hybrid_name, _one_letter
    aas = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
           "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
    for term in ("N", "C"):
        names = {hybrid_name("LEU", t, term=term) for t in aas}
        assert len(names) == len(aas), "collision among canonical terminal names"
        assert all(len(n) <= 4 for n in names), sorted(n for n in names if len(n) > 4)

