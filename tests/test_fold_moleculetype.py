"""fold_moleculetype: fold a guest moleculetype into a host so one [moleculetype]
holds both -- the covalent-complex assembly (a bond can't span moleculetypes)."""
import textwrap
from pmx.forcefield import Topology, fold_moleculetype


def _itp(name, atoms, bond, atype):
    rows = "\n".join(
        f"  {i+1}  {atype}  1  {name}  {an}  {i+1}  0.0  12.011"
        for i, an in enumerate(atoms))
    return textwrap.dedent(f"""\
        [ atomtypes ]
        ; name  at.num  mass  charge  ptype  sigma  epsilon
        {atype}  6  12.011  0.0  A  0.34  0.36
        [ moleculetype ]
        {name}  3
        [ atoms ]
        {rows}
        [ bonds ]
        {bond[0]}  {bond[1]}  1
        """)


def test_fold_two_moleculetypes(tmp_path):
    h = tmp_path / "host.itp"; h.write_text(_itp("HOS", ["A1", "A2", "A3"], (1, 2), "hh"))
    g = tmp_path / "guest.itp"; g.write_text(_itp("GST", ["B1", "B2"], (1, 2), "gg"))
    host = Topology(str(h), is_itp=True, assign_types=True, self_contained=True)
    guest = Topology(str(g), is_itp=True, assign_types=True, self_contained=True)

    na = len(host.atoms)
    fold_moleculetype(host, guest)

    # atoms concatenated; guest atoms renumbered to continue after host
    assert len(host.atoms) == 5
    assert [a.id for a in host.atoms] == [1, 2, 3, 4, 5]
    assert [a.name for a in host.atoms] == ["A1", "A2", "A3", "B1", "B2"]
    # guest residue number shifted to be disjoint from host's
    assert host.atoms[3].resnr != host.atoms[0].resnr
    # atom types merged (both hh and gg present, de-duplicated)
    names = {t["name"] for t in host.atomtypes}
    assert {"hh", "gg"} <= names
    # the guest's bond now references the renumbered atoms (ids 4,5)
    guest_bond = host.bonds[-1]
    assert {guest_bond[0].id, guest_bond[1].id} == {na + 1, na + 2}
    # two residues in the single moleculetype
    assert {r.resname for r in host.residues} == {"HOS", "GST"}
