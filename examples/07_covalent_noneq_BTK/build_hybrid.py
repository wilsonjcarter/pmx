#!/usr/bin/env python
"""Build a dual-state (A/B) GROMACS topology for an alchemical covalent-bond
formation, using pmx.covalent.edit.induced_bond_morph to define exactly which
bonded terms morph (k=0 in state A -> full GAFF value in state B).

    python build_hybrid.py  in.top  out.top  ATOM_A  ATOM_B  [heavy_only]

ATOM_A/ATOM_B are the names of the two atoms whose bond forms (e.g. SG CZ).
Pass a 5th arg ``heavy_only`` to drop induced angles/dihedrals that involve a
hydrogen -- reproducing the pre-graph disulfide hand-list, so the two topologies
differ by exactly the omitted H-torsions (used to quantify their dG effect).
Requires pmx on PYTHONPATH.  The [ atoms ]/charges are left unperturbed, so
dH/dlambda is purely the graph-derived inter-residue bonded morph.
"""
import sys, os

TOP, OUT, NA, NB = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
HEAVY_ONLY = len(sys.argv) > 5 and sys.argv[5] == "heavy_only"

sys.path.insert(0, os.environ.get("PMX_SRC",
                os.path.join(os.environ["HOME"], "Documents/pmx/src")))
from pmx.atom import Atom
from pmx.covalent.edit import induced_bond_morph

lines = open(TOP).read().splitlines()


def sections(lines):
    sec, cur = {}, None
    for ln in lines:
        s = ln.strip()
        if s.startswith("[") and s.endswith("]"):
            cur = s[1:-1].strip(); sec.setdefault(cur, [])
        elif cur is not None:
            sec[cur].append(ln)
    return sec


sec = sections(lines)


def rows(name):
    out = []
    for ln in sec.get(name, []):
        s = ln.split(";")[0].strip()
        if s:
            out.append(s.split())
    return out


names = [r[4] for r in rows("atoms")]                 # col 5 = atom name
atoms = [Atom(name=n, id=i + 1) for i, n in enumerate(names)]
for a in atoms:
    a.type = a.typeB = "X"
adj = {}
for r in rows("bonds"):
    i, j = int(r[0]), int(r[1])
    adj.setdefault(atoms[i - 1], set()).add(atoms[j - 1])
    adj.setdefault(atoms[j - 1], set()).add(atoms[i - 1])

ai = names.index(NA) + 1
bi = names.index(NB) + 1
morph = induced_bond_morph([(atoms[ai - 1], atoms[bi - 1])], adj, "form")

cross_bond = {frozenset((ai, bi))}
def ang_key(a, b, c): return (frozenset((a, c)), b)
def dih_key(q): return min(tuple(q), tuple(reversed(q)))
def _has_h(it):
    return any(a.name[:1] in "Hh" for a in it.atoms)

ang_it = [it for it in morph.interactions if it.kind == "angle"
          and not (HEAVY_ONLY and _has_h(it))]
dih_it = [it for it in morph.interactions if it.kind == "dihedral"
          and not (HEAVY_ONLY and _has_h(it))]
cross_ang = {ang_key(it.atoms[0].id, it.atoms[1].id, it.atoms[2].id) for it in ang_it}
cross_dih = {dih_key([a.id for a in it.atoms]) for it in dih_it}
print(f"forming bond: {NA}-{NB}  ->  1 bond, {len(cross_ang)} angles, "
      f"{len(cross_dih)} dihedrals morph A(k=0)->B(full)"
      + ("   [heavy-only: H-torsions dropped]" if HEAVY_ONLY else ""))


def emit_bonds():
    out = []
    for r in rows("bonds"):
        i, j = int(r[0]), int(r[1])
        if frozenset((i, j)) in cross_bond:
            out.append(f"  {i:5d} {j:5d}  {r[2]}  {r[3]} 0.0 {r[3]} {r[4]}   ; A-off/B-on")
        else:
            out.append("  " + " ".join(r))
    return out


def emit_angles():
    out = []
    for r in rows("angles"):
        i, j, k = int(r[0]), int(r[1]), int(r[2])
        if ang_key(i, j, k) in cross_ang:
            out.append(f"  {i:5d} {j:5d} {k:5d}  {r[3]}  {r[4]} 0.0 {r[4]} {r[5]}   ; A-off/B-on")
        else:
            out.append("  " + " ".join(r))
    return out


def emit_dihedrals():
    out = []
    for r in rows("dihedrals"):
        i, j, k, l = int(r[0]), int(r[1]), int(r[2]), int(r[3])
        phi, kd, mult = r[5], r[6], r[7]
        if dih_key([i, j, k, l]) in cross_dih:
            out.append(f"  {i:5d} {j:5d} {k:5d} {l:5d}  {r[4]}  {phi} 0.0 {mult} "
                       f"{phi} {kd} {mult}   ; A-off/B-on")
        else:
            out.append("  " + " ".join(r))
    return out


repl = {"bonds": emit_bonds(), "angles": emit_angles(), "dihedrals": emit_dihedrals()}
out, cur, skip = [], None, False
for ln in lines:
    s = ln.strip()
    if s.startswith("[") and s.endswith("]"):
        cur = s[1:-1].strip()
        out.append(ln)
        if cur in repl:
            out.append("; ai aj [ak al] funct  <A params>  <B params>")
            out.extend(repl[cur]); skip = True
        else:
            skip = False
        continue
    if skip:
        if s == "":
            skip = False; out.append(ln)
        continue
    out.append(ln)
open(OUT, "w").write("\n".join(out) + "\n")
print("wrote", OUT)
