#!/usr/bin/env python3
"""
Example 1 — Trp Cage W6A: Python API version (FF: charmm36m-mut)

Demonstrates the pmx Python API for:
  - loading a structure
  - applying an alchemical mutation
  - filling hybrid B states in the topology

Run: python3 run.py
(Requires mutant.pdb and topol.top from the shell workflow,
 or run after the gmx pdb2gmx step in run.sh)
"""
import os
import subprocess
from pmx.gmx import set_gmxlib
from pmx.model import Model
from pmx.forcefield import Topology
from pmx.alchemy import mutate, gen_hybrid_top

# ── 0. Environment ────────────────────────────────────────────────────────
set_gmxlib()
FF = 'charmm36m-mut'

# ── 1. Fetch structure (same as fetch_input.sh) ───────────────────────────
if not os.path.exists('input/1L2Y.pdb'):
    subprocess.run(['bash', 'input/fetch_input.sh'], check=True)

# ── 2. Load and mutate ────────────────────────────────────────────────────
print('Loading 1L2Y.pdb ...')
m = Model('input/1L2Y.pdb',
          renumber_residues=True,
          bPDBTER=True,
          rename_atoms=True,
          scale_coords='A')

print('Sequence:', ' '.join(r.resname for r in m.residues))
print('Mutating Trp6 -> Ala ...')

m2 = mutate(
    m          = m,
    mut_resid  = 6,
    mut_resname= 'A',       # one-letter code for Ala
    ff         = FF,
    inplace    = False,
    verbose    = True,
)

m2.write('mutant.pdb')
print('Written: mutant.pdb')

# ── 3. Run pdb2gmx (requires GROMACS) ────────────────────────────────────
print('\nRunning gmx pdb2gmx ...')
ret = subprocess.run([
    'gmx', 'pdb2gmx',
    '-f', 'mutant.pdb',
    '-o', 'processed.gro',
    '-p', 'topol.top',
    '-ff', FF,
    '-water', 'tip3p',
    '-ignh',
], capture_output=True, text=True)
if ret.returncode != 0:
    print('pdb2gmx stdout:', ret.stdout[-2000:])
    print('pdb2gmx stderr:', ret.stderr[-2000:])
    raise RuntimeError('gmx pdb2gmx failed')
print('Written: processed.gro, topol.top')

# ── 4. Fill hybrid B states ───────────────────────────────────────────────
print('\nFilling B states ...')
topol = Topology('topol.top', ff=FF, version='new')
pmxtop, pmxitps = gen_hybrid_top(topol=topol, recursive=True, verbose=True)
pmxtop.write('pmxtop.top')
print('Written: pmxtop.top')

# ── 5. Sanity check ───────────────────────────────────────────────────────
qA = pmxtop.get_qA()
qB = pmxtop.get_qB()
print(f'\nState A charge: {qA:.1f} e')
print(f'State B charge: {qB:.1f} e')
print(f'ΔQ (A→B):       {qB - qA:.1f} e  (should be 0 for W6A)')

# Confirm the hybrid residue is in the topology
hybrid_residues = [r for r in pmxtop.residues if r.is_hybrid()]
print(f'\nHybrid residues in pmxtop.top:')
for r in hybrid_residues:
    print(f'  Residue {r.id}: {r.resname}')
