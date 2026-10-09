"""Build a real GACE-CYS-GLY (glutathione) structure from the force-field RTP --
GACE is the N-terminal gamma-glutamyl residue (linked to CYS through its
gamma-carbonyl), NOT a modified ACE cap.  Rough peptide-bond placement is enough;
pdb2gmx forms the inter-residue bonds from the rtp."""
import sys, numpy as np
from pmx.model import Model
from pmx.chain import Chain
from pmx.covalent.target import TargetResidue

ff, out = sys.argv[1], sys.argv[2]
# all internal forms -- let pdb2gmx -ter add/name the terminal atoms per the ff
# (building GLY as term='C' bakes in an OC2 that some ffs' CGLY rtp don't use).
SEQ = [("GACE", None), ("CYS", None), ("GLY", None), ("NME", None)]

def atom(res, name):
    for a in res.atoms:
        if a.name == name:
            return a
    return None

residues = []
for i, (name, term) in enumerate(SEQ):
    r = TargetResidue.from_rtp(name, ff, term=term).residue
    r.set_resname(name)
    residues.append(r)

# link residue i's C (carbonyl) to residue i+1's N: translate i+1 so N sits ~1.33 A
# from the previous C, roughly opposite the C=O direction.
for i in range(1, len(residues)):
    prevC = atom(residues[i-1], "C")
    N = atom(residues[i], "N")
    if prevC is None or N is None:
        raise SystemExit(f"missing C/N linking {SEQ[i-1][0]}-{SEQ[i][0]}")
    prevO = atom(residues[i-1], "O")
    d = np.array(prevC.x) - np.array(prevO.x); d = d/np.linalg.norm(d)
    target = np.array(prevC.x) + 1.33*d
    shift = target - np.array(N.x)
    for a in residues[i].atoms:
        a.x = list(np.array(a.x) + shift)

ch = Chain()
ch.residues = residues
for i, r in enumerate(residues):
    r.set_resid(i+1); r.set_chain_id("A"); r.chain = ch
    for a in r.atoms:
        a.chain_id = "A"
ch.set_chain_id("A")
m = Model()
m.chains = [ch]; ch.model = m
m.al_from_resl(); m.resl_from_chains() if hasattr(m,'resl_from_chains') else None
m.residues = residues
m.al_from_resl(); m.renumber_atoms()
m.write(out)
print("wrote", out, "with", [r.resname for r in residues], "natoms", len(m.atoms))
