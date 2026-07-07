"""Ligand atomtypes supplied via a #include'd ITP (not injected into the FF).

pmx's NBParser must read EVERY [ atomtypes ] block of the (cpp-expanded) topology,
so a ligand ITP can contribute a second block on top of the force field's -- and
must ignore stray banner/blank lines that a multi-block read can capture.
"""
from pmx.parser import readSectionAll, readSection
from pmx.ffparser import NBParser


def test_readSectionAll_collects_every_block_and_skips_blanks():
    lines = [
        "[ atomtypes ]\n", "A 6 12.0 0 A 0.3 0.4\n", "\n",
        "[ bondtypes ]\n", "A A 1 0.1 100\n",
        "[ atomtypes ]\n", "B 6 12.0 0 A 0.3 0.4\n",
    ]
    # readSection stops at the first block; readSectionAll spans both
    first = [l.split()[0] for l in readSection(lines, "[ atomtypes ]", "[") if l.strip()]
    assert first == ["A"]
    got = [l.split()[0] for l in readSectionAll(lines, "[ atomtypes ]", "[") if l.strip()]
    assert got == ["A", "B"]


def test_nbparser_reads_ff_plus_ligand_blocks_ignoring_banner():
    # a force-field .itp style banner (raw '****', not a ';' comment) + two
    # [ atomtypes ] blocks: the FF's (incl. an amber digit-leading type '2C')
    # and a ligand ITP's GAFF block.
    fulltop = [
        "********************************************************************\n",
        "* banner text *\n",
        "********************************************************************\n",
        "[ atomtypes ]\n",
        "CX 6 12.01 0.0 A 3.39967e-01 3.59824e-01\n",
        "2C 6 12.01 0.0 A 3.39967e-01 3.59824e-01\n",
        "\n",
        "[ bondtypes ]\n",
        "CX 2C 1 0.15 200000\n",
        "; ligand ITP contributes a second block\n",
        "[ atomtypes ]\n",
        "c3 6 12.010 0.0 A 0.33977 0.45103\n",
        "ca 6 12.010 0.0 A 0.33152 0.41338\n",
    ]
    nb = NBParser(fulltop, version="new", ff="amber14sbmut")
    # both the FF's and the ligand's types are known...
    for t in ("CX", "2C", "c3", "ca"):
        assert t in nb.atomtypes, "%s missing from parsed atomtypes" % t
    # ...the digit-leading amber type survived the data-line filter...
    assert nb.atomtypes["2C"]["bond_type"] == "2C"
    # ...and the '****' banner did not become a bogus atomtype.
    assert not any(k.startswith("*") for k in nb.atomtypes)
