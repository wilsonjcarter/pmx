#!/usr/bin/env python3
"""
Example 5 — Phosphoserine (SEP): Python API version

Full NSAA workflow:
  1. prepare_nsaa_ff  ->  SEP.rtp, SEP_bonded.itp
  2. mutate_nsaa      ->  mutant.pdb, S2SEP.mtp
  3. pdb2gmx          ->  topol.top
  4. gen_hybrid_top   ->  pmxtop.top

Run: python3 run.py
(Requires input/protein.pdb, input/sep.pdb, input/sep.itp)
"""
import os
import subprocess
from pmx.gmx import set_gmxlib
from pmx.model import Model
from pmx.forcefield import Topology
from pmx.alchemy import gen_hybrid_top
from pmx.nsaa import mutate_nsaa, prepare_nsaa_ff
from pmx.utils import get_ff_path

# ── Configuration ─────────────────────────────────────────────────────────
FF         = 'charmm36m-mut'
WATER      = 'tip3p'
MUT_RESID  = 42       # Ser residue to phosphorylate (pmx-renumbered, 1-based)
MUT_RESNAME= 'SEP'   # target NSAA name
BASE_RES   = 'SER'   # canonical AA the NSAA is derived from

# ── 0. Environment ────────────────────────────────────────────────────────
set_gmxlib()
ffdir = get_ff_path(FF)
print(f'Force field: {FF}')
print(f'FF directory: {ffdir}')

# ── 1. prepare_nsaa_ff ────────────────────────────────────────────────────
# Reads sep.itp (GAFF2-parameterised SEP fragment with CH3 cap at CA),
# builds a proper RTP entry for SEP, and writes the bonded supplement ITP.
print('\n[Step 1] prepare_nsaa_ff ...')
out_rtp, out_supplement = prepare_nsaa_ff(
    itp_file       = 'input/sep.itp',
    base_resname   = BASE_RES,
    new_resname    = MUT_RESNAME,
    ff             = FF,
    target_charge  = -2,      # net charge of phosphoserine at pH 7
    atomtypes_itp  = None,    # set if your ITP uses abstract types
    out_rtp        = 'SEP.rtp',
    out_supplement = 'SEP_bonded.itp',
)
print(f'  RTP patch       : {out_rtp}')
print(f'  Bonded supplement: {out_supplement}')

# ── 2. mutate_nsaa ────────────────────────────────────────────────────────
# Loads the protein, builds the hybrid structure (Ser sidechain in state A,
# SEP sidechain in state B), writes mutant.pdb and S2SEP.mtp.
print('\n[Step 2] mutate_nsaa ...')
m = Model('input/protein.pdb',
          renumber_residues=True,
          bPDBTER=True,
          rename_atoms=True,
          scale_coords='A')

print(f'  Loaded {len(m.residues)} residues, {len(m.atoms)} atoms')
print(f'  Mutating residue {MUT_RESID} ({m.residues[MUT_RESID-1].resname}) -> {MUT_RESNAME}')

mutate_nsaa(
    m            = m,
    mut_resid    = MUT_RESID,
    mut_resname  = MUT_RESNAME,
    ff           = FF,
    itp          = 'input/sep.itp',
    nonstandardPDB = 'input/sep.pdb',
    inplace      = True,
    verbose      = True,
)
m.write('mutant.pdb')
print('  Written: mutant.pdb')

# Expected MTP output name: S2SEP.mtp (src one-letter + '2' + resname)
mtp_file = f'S2{MUT_RESNAME}.mtp'
if not os.path.exists(mtp_file):
    raise FileNotFoundError(
        f'Expected MTP file {mtp_file} not found. '
        'Check that mutate_nsaa wrote it to the working directory.')
print(f'  MTP file: {mtp_file}')

# ── 3. pdb2gmx ────────────────────────────────────────────────────────────
# pdb2gmx uses the SEP.rtp entry written into the FF directory by Step 1.
print('\n[Step 3] gmx pdb2gmx ...')
ret = subprocess.run([
    'gmx', 'pdb2gmx',
    '-f', 'mutant.pdb',
    '-o', 'processed.gro',
    '-p', 'topol.top',
    '-ff', FF,
    '-water', WATER,
    '-ignh',
], capture_output=True, text=True)
if ret.returncode != 0:
    print('pdb2gmx stderr:', ret.stderr[-3000:])
    raise RuntimeError('gmx pdb2gmx failed')
print('  Written: processed.gro, topol.top')

# ── 4. gen_hybrid_top ────────────────────────────────────────────────────
# extra_mtp_files  : tells gentop about the S2SEP hybrid residue
# supplement_bonded_files: injects GAFF/cross-boundary bonded terms before
#                          the B-state parameter lookup
print('\n[Step 4] gen_hybrid_top ...')
topol = Topology('topol.top', ff=FF, version='new')
pmxtop, pmxitps = gen_hybrid_top(
    topol                  = topol,
    recursive              = True,
    verbose                = True,
    extra_mtp_files        = [mtp_file],
    supplement_bonded_files= [out_supplement],
)
pmxtop.write('pmxtop.top')
print('  Written: pmxtop.top')

# ── 5. Summary ───────────────────────────────────────────────────────────
qA = pmxtop.get_qA()
qB = pmxtop.get_qB()
hybrid_res = [r for r in pmxtop.residues if r.is_hybrid()]
print(f'\n=== Done ===')
print(f'State A charge : {qA:.0f} e')
print(f'State B charge : {qB:.0f} e  (phosphate: delta = {qB-qA:.0f} e)')
print(f'Hybrid residues: {[f"{r.id}:{r.resname}" for r in hybrid_res]}')
print()
print(f'NSAA FF patch  : {ffdir}/SEP.rtp')
print(f'Bonded suppl.  : SEP_bonded.itp')
print(f'Hybrid MTP     : {mtp_file}')
print(f'Hybrid struct  : mutant.pdb')
print(f'Hybrid topology: pmxtop.top')
print()
print('Next: solvate, add 2 Na+ counter-ions for state B, em/eq, production.')
