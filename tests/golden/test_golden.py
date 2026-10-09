"""Golden-output guard for the mutation/hybrid engine.

Regenerates the hybrid structure + dual-state topology for a representative matrix
of cases (point mutations, protonations, terminal edit) and requires them
byte-identical to the committed fixtures.  This is the safety net for the
mutation-engine cleanup (see MUTATION_ENGINE_CLEANUP.md): any refactor that changes
an emitted hybrid fails here.

Requires ``gmx`` on PATH and ``GMXLIB`` pointing at a force-field directory that
contains ``charmm36m-mut.ff`` (the canonical ``pmx/data/mutff`` is used by default
if GMXLIB is unset).
"""
import os
import shutil
import subprocess
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
GEN = os.path.join(HERE, "generate_golden.py")

pytestmark = pytest.mark.skipif(
    shutil.which("gmx") is None, reason="gmx not on PATH")


def test_golden_hybrids_unchanged():
    env = dict(os.environ)
    env.setdefault("GMXLIB", os.path.join(HERE, "..", "..", "src", "pmx", "data", "mutff"))
    env["PYTHONPATH"] = os.path.join(HERE, "..", "..", "src") + os.pathsep + env.get("PYTHONPATH", "")
    r = subprocess.run([sys.executable, GEN, "--check"], env=env,
                       capture_output=True, text=True)
    drift = [l for l in r.stdout.splitlines() if l.startswith("DRIFT")]
    assert not drift and r.returncode == 0, (
        "golden hybrid output drifted:\n" + "\n".join(drift) + "\n" + r.stdout)
