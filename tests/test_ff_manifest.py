#!/usr/bin/env python
"""Guard the hand-introduced force-field content declared in pmx_ff.manifest.

Every residue and parameter block pmx adds to a bundled force field lives inside
files that are also vendored or regenerated.  On 2026-07-15 a revert of
``charmm36m-mut/merged.rtp`` silently removed the CN0/KN0/KS0 NOS endpoints, 7
acyl-lysines and NME; nothing failed until weeks later, and recovery meant
digging the parameters back out of build artifacts.

These tests turn that class of loss into a failing test:

* :func:`test_manifest_entries_resolve` -- everything declared still resolves.
* :func:`test_no_regression_against_upstream` -- no *base RTP* entry present
  upstream has gone missing here (skipped when the upstream checkout is absent).

Run standalone for a readable report::

    python tests/test_ff_manifest.py
"""
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
MUTFF = os.path.join(HERE, "..", "src", "pmx", "data", "mutff")
MANIFEST = os.path.join(MUTFF, "pmx_ff.manifest")

# Read-only reference checkout; absent on other machines, so those checks skip.
UPSTREAM = os.path.expanduser(
    "~/Documents/GitHub/pmx/src/pmx/data/ff/mutff")


def _records():
    """Parse the manifest into (ff, kind, value, parent) tuples."""
    out = []
    with open(MANIFEST) as fh:
        for line in fh:
            line = line.split("#")[0].strip()
            if not line:
                continue
            f = line.split()
            if len(f) < 3:
                continue
            out.append((f[0], f[1], f[2], f[3] if len(f) > 3 else None))
    return out


def _ffdir(ff):
    return os.path.join(MUTFF, ff + ".ff")


def _check_all():
    """Return a list of human-readable failures ([] when healthy)."""
    os.environ.setdefault("GMXLIB", os.path.abspath(MUTFF))
    from pmx.ffparser import load_ff_rtp

    fails = []
    rtp_cache = {}
    for ff, kind, value, parent in _records():
        if ff != "-" and not os.path.isdir(_ffdir(ff)):
            fails.append("%s: force-field directory missing" % ff)
            continue

        if kind == "residue":
            try:
                if ff not in rtp_cache:
                    rtp_cache[ff] = load_ff_rtp(ff)
                rtp = rtp_cache[ff]
            except Exception as e:                      # pragma: no cover
                fails.append("%s: could not load RTP database (%s)" % (ff, e))
                continue
            if value not in rtp:
                searched = ", ".join(os.path.basename(s) for s in rtp.sources)
                fails.append(
                    "%s: residue %s MISSING (expected in %s; searched %s)"
                    % (ff, value, parent or "?", searched))

        elif kind == "charge":
            # `parent` holds the expected net charge.  Guards the ECC scaling:
            # charmm36m-ecc-mut shipped SEP/TPO/PTR/MLY/M3L at full formal charge
            # while its canonical residues were scaled, so any ECC phospho or
            # methyl calculation mixed scaled protein charges with an unscaled
            # modification.  Scaled 2026-08-06; this stops it drifting back.
            try:
                if ff not in rtp_cache:
                    rtp_cache[ff] = load_ff_rtp(ff)
                rtp = rtp_cache[ff]
            except Exception as e:                      # pragma: no cover
                fails.append("%s: could not load RTP database (%s)" % (ff, e))
                continue
            if value not in rtp:
                fails.append("%s: residue %s MISSING (charge check)" % (ff, value))
                continue
            q = sum(float(a[2]) for a in rtp[value]["atoms"])
            want = float(parent)
            if abs(q - want) > 1e-4:
                fails.append("%s: %s net charge %+.4f, expected %+.4f"
                             % (ff, value, q, want))

        elif kind == "hdb":
            found = False
            for name in sorted(os.listdir(_ffdir(ff))):
                if not name.endswith(".hdb"):
                    continue
                with open(os.path.join(_ffdir(ff), name), errors="ignore") as fh:
                    for line in fh:
                        w = line.split()
                        if w and w[0] == value:
                            found = True
                            break
                if found:
                    break
            if not found:
                fails.append("%s: no hydrogen-build (.hdb) entry for %s" % (ff, value))

        elif kind == "file":
            if not os.path.isfile(os.path.join(_ffdir(ff), value)):
                fails.append("%s: file %s MISSING" % (ff, value))

        elif kind == "include":
            path = os.path.join(_ffdir(ff), value)
            if not os.path.isfile(path):
                fails.append("%s: include target %s MISSING" % (ff, value))
                continue
            parent_path = os.path.join(_ffdir(ff), parent or "")
            if not os.path.isfile(parent_path):
                fails.append("%s: include parent %s MISSING" % (ff, parent))
                continue
            with open(parent_path, errors="ignore") as fh:
                if value not in fh.read():
                    fails.append("%s: %s exists but is not #include'd by %s"
                                 % (ff, value, parent))

        elif kind == "restype":
            # A per-force-field residuetypes.dat REPLACES the tree-level one.
            strays = [d for d in sorted(os.listdir(MUTFF))
                      if d.endswith(".ff")
                      and os.path.isfile(os.path.join(MUTFF, d, "residuetypes.dat"))]
            if strays:
                fails.append(
                    "per-force-field residuetypes.dat shadows the tree-level one in: %s"
                    % ", ".join(strays))
            root = os.path.join(MUTFF, "residuetypes.dat")
            with open(root, errors="ignore") as fh:
                names = {l.split()[0] for l in fh if l.split()}
            if value not in names:
                fails.append("residuetypes.dat: %s not classified" % value)

    # de-duplicate (the residuetypes stray check repeats per record)
    seen, uniq = set(), []
    for f in fails:
        if f not in seen:
            seen.add(f)
            uniq.append(f)
    return uniq


def test_manifest_entries_resolve():
    fails = _check_all()
    assert not fails, ("declared force-field content is missing:\n  "
                       + "\n  ".join(fails))


def _all_entries(path):
    """Union of residue names across a force-field directory's non-hybrid RTPs.

    Deliberately not just the base RTP: an entry legitimately *moves* from the
    vendored file into a topic file (acyl.rtp, nos.rtp, ...), which is the whole
    point of the split and must not read as a loss.  Hybrid libraries are
    excluded, matching ffparser.load_ff_rtp.
    """
    from pmx.ffparser import RTPParser, _is_hybrid_rtp
    names = set()
    for f in sorted(os.listdir(path)):
        if not f.endswith(".rtp") or _is_hybrid_rtp(f):
            continue
        try:
            names |= set(RTPParser(os.path.join(path, f)).keys)
        except Exception:                                # pragma: no cover
            continue
    return names


@pytest.mark.skipif(not os.path.isdir(UPSTREAM),
                    reason="upstream reference checkout not present")
def test_no_regression_against_upstream():
    """No residue that upstream ships may go missing here.

    This is the check that would have caught the 2026-07-15 loss the day it
    happened: 8 entries silently vanished from charmm36m-mut/merged.rtp.
    Entries we *add* are fine -- only disappearance is a failure.
    """
    lost = {}
    for ff in sorted(d for d in os.listdir(MUTFF) if d.endswith(".ff")):
        up = os.path.join(UPSTREAM, ff)
        if not os.path.isdir(up):
            continue
        try:
            missing = _all_entries(up) - _all_entries(os.path.join(MUTFF, ff))
        except Exception:                                # pragma: no cover
            continue
        if missing:
            lost[ff] = sorted(missing)
    assert not lost, ("residues present upstream but missing here:\n  "
                      + "\n  ".join("%s: %s" % (k, ", ".join(v))
                                    for k, v in lost.items()))


if __name__ == "__main__":
    sys.path.insert(0, os.path.join(HERE, "..", "src"))
    problems = _check_all()
    print("pmx ff doctor -- %d manifest records" % len(_records()))
    if problems:
        print("\nPROBLEMS (%d):" % len(problems))
        for p in problems:
            print("  " + p)
    else:
        print("  all declared content present")
    if os.path.isdir(UPSTREAM):
        lost = {}
        for ff in sorted(d for d in os.listdir(MUTFF) if d.endswith(".ff")):
            up = os.path.join(UPSTREAM, ff)
            if os.path.isdir(up):
                try:
                    m = _all_entries(up) - _all_entries(os.path.join(MUTFF, ff))
                except Exception:
                    continue
                if m:
                    lost[ff] = sorted(m)
        print("\nupstream drift:")
        if lost:
            for k, v in lost.items():
                print("  %s: LOST %s" % (k, ", ".join(v)))
        else:
            print("  no residues lost relative to upstream")
    else:
        print("\nupstream drift: skipped (reference checkout not present)")
    sys.exit(1 if problems else 0)
