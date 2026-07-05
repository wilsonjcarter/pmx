#!/usr/bin/env python
"""TargetResidue.from_itp: build a target residue (with geometry) from an ITP
fragment alone -- RDKit-embedded coordinates from the bond graph, like from_rtp,
so no separate structure file is needed."""
import textwrap
import pytest

pytest.importorskip("rdkit")   # from_itp embeds coordinates via RDKit

from pmx.covalent import TargetResidue

_FRAG_ITP = textwrap.dedent("""\
    [ atoms ]
    ;  nr  type  resnr res  name cgnr  charge   mass
       1   CT3   1     FRG  CA   1     0.000   12.011
       2   HA3   1     FRG  HA1  1     0.090    1.008
       3   HA3   1     FRG  HA2  1     0.090    1.008
       4   HA3   1     FRG  HA3  1     0.090    1.008
       5   CT2   1     FRG  CB   2     0.050   12.011
       6   HA2   1     FRG  HB1  2     0.090    1.008
       7   HA2   1     FRG  HB2  2     0.090    1.008
       8   OH1   1     FRG  OG   3    -0.660   15.999
       9   H     1     FRG  HG   3     0.430    1.008
    [ bonds ]
       1  2
       1  3
       1  4
       1  5
       5  6
       5  7
       5  8
       8  9
""")


def test_from_itp_builds_residue_with_embedded_geometry(tmp_path):
    itp = tmp_path / "frag.itp"
    itp.write_text(_FRAG_ITP)

    res = TargetResidue.from_itp(str(itp), "SFR").residue

    # atoms, types, charges taken from the ITP
    assert [a.name for a in res.atoms] == \
        ["CA", "HA1", "HA2", "HA3", "CB", "HB1", "HB2", "OG", "HG"]
    assert res.atoms[7].atomtype == "OH1"           # OG
    assert abs(res.atoms[7].q - (-0.660)) < 1e-6

    # geometry embedded from the bond graph (no structure file) -> real 3-D coords
    assert any(any(abs(c) > 1e-6 for c in a.x) for a in res.atoms)

    # connectivity attached from the ITP [ bonds ]
    cb = next(a for a in res.atoms if a.name == "CB")
    assert {b.name for b in cb.bonds} == {"CA", "HB1", "HB2", "OG"}
