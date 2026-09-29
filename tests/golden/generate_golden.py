#!/usr/bin/env python
"""Generate the golden hybrid-topology fixtures.

For each case in MATRIX, mutate the shared test peptide (``_peptide.pdb``) and
capture, normalized for stable diffing:
  * ``mutant.pdb``  -- the hybrid structure the engine emits
  * ``hybrid.top``  -- the dual-state topology after pdb2gmx + gentop (the artifact
                       free-energy actually consumes)

Run with ``--check`` to regenerate into a temp dir and diff against the committed
fixtures instead of overwriting them (this is what the pytest does).

Determinism: mutate/gentop are deterministic; RDKit embedding is seeded. gmx is run
with a fixed FF (GMXLIB) so includes are stable; normalization strips paths/comments.

Usage:
    GMXLIB=<ff dir> python generate_golden.py            # (re)write fixtures
    GMXLIB=<ff dir> python generate_golden.py --check     # diff-only, exit 1 on drift
"""
import os, sys, shutil, subprocess, tempfile, re, glob
from copy import deepcopy

# A runtime-built hybrid residue is named <one_letter>2<...> (e.g. K2AL) and the
# generative engine accumulates hybrids into nsaa.{rtp,mtp,arn}.  Both are engine
# *output*, not part of the committed force field -- so a working tree can carry
# leftovers from earlier runs.  The isolated generative fixtures must start from a
# clean FF (else pdb2gmx sees the same residue in a stale <hybrid>.rtp AND the
# freshly written nsaa.rtp -> "defined in two places").  No shipped FF file has a
# 4-char, '2'-in-second-position stem, so this pattern only ever matches hybrids.
_HYBRID_STEM = re.compile(r'^[A-Za-z]2[A-Za-z]{2}$')


def _strip_runtime_hybrids(ffdir):
    """Remove engine-generated hybrid artifacts from a *copied* FF directory."""
    for f in os.listdir(ffdir):
        stem, ext = os.path.splitext(f)
        if ext in (".rtp", ".mtp", ".arn") and (stem == "nsaa"
                                                 or _HYBRID_STEM.match(stem)):
            os.remove(os.path.join(ffdir, f))

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
from pmx.model import Model
from pmx.alchemy import mutate
from pmx.nsaa import mutate_nsaa
from pmx.covalent import amidate_cterminus

PEPTIDE = os.path.join(HERE, "_peptide.pdb")
# a second input whose cysteine is INTERNAL (PHE ARG THR LEU LYS ASN CYS TRP GLN)
CYS_PEPTIDE = os.path.join(HERE, "..", "data", "chain_create_FRTLKNCWQ.pdb")
FF = "charmm36m-mut"
ALY_TMPL = os.path.join(HERE, "aly_lig.pdb")   # a complete acetyl-Lys template
SFR_FRAG = os.path.join(HERE, "sfr_frag.itp")  # a CH3-capped side-chain ITP fragment

# case name -> (kind, kwargs).  resid is the peptide position (1..9).
#   ALA1 VAL2 ASP3 GLU4 LEU5 ILE6 LYS7 HSD8 CYS9
MATRIX = [
    # -- point mutations --
    ("point__ala2val",   ("mutate", dict(resid=1, target="VAL"))),
    ("point__val2ala",   ("mutate", dict(resid=2, target="ALA"))),
    ("point__asp2lys",   ("mutate", dict(resid=3, target="LYS"))),
    ("point__leu2phe",   ("mutate", dict(resid=5, target="PHE"))),   # ring symmetry
    ("point__ile2val",   ("mutate", dict(resid=6, target="VAL"))),   # methyl symmetry
    # one-letter target code at a terminus -> generative path (regression guard for
    # the one-letter resolution; three-letter equivalents are the cases above).
    ("point__cterm_1letter", ("mutate", dict(resid=9, target="A"))),  # CYS9->ALA, C-term
    # -- protonations (charge-changing; same machinery) --
    ("proton__asp2aspp", ("mutate", dict(resid=3, target="ASPP"))),
    ("proton__glu2glup", ("mutate", dict(resid=4, target="GLUP"))),
    ("proton__hsd2hsp",  ("mutate", dict(resid=8, target="HSP"))),
    ("proton__lys2lyn",  ("mutate", dict(resid=7, target="LYN"))),
    # CYS->CYM on an INTERNAL cysteine.  The shared peptide carries CYS only at
    # position 9 (its C-terminus), and a terminal CYS->CYM would need the 5-char
    # hybrid name C2CMc, which cannot round-trip through a PDB and is refused --
    # so this case uses a peptide with an internal cysteine (FRTLKNCWQ, CYS7).
    ("proton__cys2cym",  ("mutate", dict(resid=7, target="CYM",
                                        peptide=CYS_PEPTIDE))),
    # -- NSAA (Stage-1 sources; these WRITE the hybrid into the FF, so build_case
    #    isolates them in a temp FF copy to avoid polluting the tracked force field) --
    ("nsaa__lys2aly_struct", ("mutate_nsaa", dict(resid=7, resname="ALY", nsapdb=ALY_TMPL))),
    ("nsaa__lys2aly_rtp",    ("mutate_gen",  dict(resid=7, target="ALY"))),
    # from_library front end: a standard AA built by build_chain (nsapdb=None) and
    # mapped by the nsaa engine (domapping) -- guards mutate_nsaa's build_chain path.
    ("nsaa__leu2phe_lib",    ("mutate_nsaa", dict(resid=5, resname="PHE", nsapdb=None))),
    # from_itp front end: a CH3-capped fragment (no structure) -- guards the
    # capped-ITP path.  `pmx mutate_nsaa --itp` now delegates to the covalent
    # engine (geometry embedded from the ITP graph), so this also guards that
    # delegation producing the same hybrid as `pmx mutate --itp`.
    ("nsaa__leu2sfr_itp",    ("mutate_nsaa", dict(resid=5, resname="SFR", nsapdb=None, itp=SFR_FRAG))),
    # -- terminal --
    ("term__cterm_amide", ("amidate", dict(resid=9))),
    ("term__cterm_del",   ("delete", dict(resid=9))),   # delete CYS9 (neighbour HSD8 -> C-term)
    ("term__nterm_del",   ("delete", dict(resid=1))),   # delete ALA1 (neighbour VAL2 -> N-term)
]


def norm_top(text):
    """Strip machine-specific lines so the topology diffs stably."""
    out = []
    for l in text.splitlines():
        s = l.strip()
        if s.startswith(";"):                       # comments / banners
            continue
        if "#include" in l and ("/" in l):          # FF include paths
            out.append(re.sub(r'"[^"]*/([^/"]+)"', r'"\1"', l))
            continue
        if not s:
            continue
        out.append(l.rstrip())
    return "\n".join(out) + "\n"


def norm_pdb(text):
    """Keep only ATOM/HETATM records with fields that matter (name/res/resid/xyz)."""
    out = []
    for l in text.splitlines():
        if l[:6] in ("ATOM  ", "HETATM"):
            out.append(l[12:26] + l[30:54])         # name+res+chain+resid + xyz
    return "\n".join(out) + "\n"


def build_case(name, kind, kwargs, outdir):
    os.makedirs(outdir, exist_ok=True)
    work = tempfile.mkdtemp()
    cwd = os.getcwd()
    saved_gmxlib = os.environ.get("GMXLIB")
    try:
        os.chdir(work)
        # EVERY case gets its own throwaway force-field directory.
        #
        # `mutate` is the generative facade and `amidate` accumulates too, so any
        # case can write hybrids into <GMXLIB>/<ff>.ff.  Isolating only some of
        # them left the harness writing into the TRACKED force field and made it
        # non-idempotent: the cases run sequentially, so each one built against a
        # force field the earlier ones had mutated.  Measured on a clean tree,
        # run 1 gave ok=14/no-top=4 and run 2 ok=4/no-top=14, and which case
        # drifted moved between runs.  A fixture generated from a clean force
        # field could therefore never match a --check run.
        #
        # This is a TEST-HARNESS measure only.  Ordinary use is unaffected: a
        # user who wants hybrids kept beside their project passes `-ff <dir>`,
        # which get_ff_path() accepts (it checks os.path.isdir first), and the
        # engine then writes there instead of into the installation.
        src_ff = os.path.join(saved_gmxlib or "", FF + ".ff")
        tmp_gmxlib = os.path.join(work, "ff")
        os.makedirs(tmp_gmxlib)
        shutil.copytree(src_ff, os.path.join(tmp_gmxlib, FF + ".ff"))
        _strip_runtime_hybrids(os.path.join(tmp_gmxlib, FF + ".ff"))
        # mirror the whole top-level GMXLIB: copying only residuetypes.dat left the
        # isolated run unable to resolve the water model, so pdb2gmx omitted the
        # `#include "tip3p.itp"` line and every fixture showed a one-line drift.
        for extra in sorted(os.listdir(saved_gmxlib or ".")):
            q = os.path.join(saved_gmxlib or "", extra)
            if os.path.isfile(q):
                shutil.copy(q, tmp_gmxlib)
        os.environ["GMXLIB"] = tmp_gmxlib
        m = Model(kwargs.get("peptide", PEPTIDE), rename_atoms=True)
        if kind == "mutate":
            mut = mutate(m, mut_resid=kwargs["resid"], mut_resname=kwargs["target"], ff=FF)
        elif kind == "mutate_gen":
            mut = mutate(m, mut_resid=kwargs["resid"], mut_resname=kwargs["target"],
                         ff=FF, generative=True)
        elif kind == "mutate_nsaa":
            mutate_nsaa(m, mut_resid=kwargs["resid"], mut_resname=kwargs["resname"],
                        ff=FF, nonstandardPDB=kwargs["nsapdb"],
                        itp=kwargs.get("itp"), inplace=True)
            mut = m
        elif kind == "amidate":
            mut = amidate_cterminus(m, resid=kwargs["resid"], ff=FF)
        elif kind == "delete":
            # terminal deletion: the neighbour morphs internal->terminus (XdeC/XdeN
            # hybrid, accumulated into nsaa.*), the target dummifies at gentop.
            mut = mutate(m, mut_resid=kwargs["resid"], mut_resname="DEL", ff=FF)
        else:
            raise ValueError(kind)
        mut.write("mutant.pdb")
        # pdb2gmx + gentop -> dual-state topology
        subprocess.run(["gmx", "pdb2gmx", "-f", "mutant.pdb", "-ff", FF, "-water",
                        "tip3p", "-o", "hyb.pdb", "-p", "hyb.top"],
                       stdin=subprocess.DEVNULL, capture_output=True)
        # nsaa.mutate_nsaa writes the hybrid MTP as <hybrid>.mtp in the working
        # dir (its documented usage is to pass it to gen_hybrid_top); covalent.mutate
        # instead accumulates into the FF's nsaa.mtp, which gentop auto-discovers.
        # Feed any working-dir MTP so both emit paths resolve.
        gentop_cmd = ["pmx", "gentop", "-p", "hyb.top", "-o", "fe.top", "-ff", FF]
        for extra_mtp in sorted(glob.glob("*.mtp")):
            gentop_cmd += ["--extra_mtp", extra_mtp]
        gt = subprocess.run(gentop_cmd, capture_output=True, text=True)
        pdb = norm_pdb(open("mutant.pdb").read())
        open(os.path.join(outdir, "mutant.pdb"), "w").write(pdb)
        if os.path.isfile("fe.top"):
            open(os.path.join(outdir, "hybrid.top"), "w").write(norm_top(open("fe.top").read()))
            return "ok"
        return "no-top (%s)" % (gt.stderr.strip().splitlines()[-1] if gt.stderr else "?")
    except Exception as e:
        return "ERROR %s" % e
    finally:
        os.chdir(cwd)
        if saved_gmxlib is not None:
            os.environ["GMXLIB"] = saved_gmxlib
        shutil.rmtree(work, ignore_errors=True)


def main():
    check = "--check" in sys.argv
    root = tempfile.mkdtemp() if check else HERE
    drift = 0
    for name, (kind, kwargs) in MATRIX:
        outdir = os.path.join(root, name)
        status = build_case(name, kind, kwargs, outdir)
        if check:
            gold = os.path.join(HERE, name)
            for f in ("mutant.pdb", "hybrid.top"):
                g, n = os.path.join(gold, f), os.path.join(outdir, f)
                if not os.path.isfile(g):
                    continue                       # no committed fixture to guard
                if not os.path.isfile(n):
                    # A committed fixture with no regenerated counterpart means
                    # the engine stopped producing that artifact -- a regression,
                    # not a pass.  (Guards the case where a refactor breaks the
                    # emit/gentop step: silently skipping it would hide the break.)
                    print("DRIFT %s/%s (not regenerated)" % (name, f)); drift += 1
                elif open(g).read() != open(n).read():
                    print("DRIFT %s/%s" % (name, f)); drift += 1
        print("%-24s %s" % (name, status))
    if check:
        shutil.rmtree(root, ignore_errors=True)
        sys.exit(1 if drift else 0)


if __name__ == "__main__":
    main()
