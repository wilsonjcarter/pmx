#!/usr/bin/env python

# pmx  Copyright Notice
# ============================
#
# The pmx source code is copyrighted, but you can freely use and
# copy it as long as you don't change or remove any of the copyright
# notices.
#
# ----------------------------------------------------------------------
# pmx is Copyright (C) 2006-2013 by Daniel Seeliger
#
#                        All Rights Reserved
#
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name of Daniel Seeliger not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
#
# DANIEL SEELIGER DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
# SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS.  IN NO EVENT SHALL DANIEL SEELIGER BE LIABLE FOR ANY
# SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
# RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
# CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
# ----------------------------------------------------------------------

"""Build hybrid structure files for non-standard amino acid (NSAA / PTM)
free energy calculations.

Two-step workflow
-----------------
**Step 1** — build the hybrid structure and write FF input files::

    pmx mutate_nsaa -f protein.pdb -o mutant.pdb \\
        --resid 42 --resname SEP --ff charmm36m-mut \\
        --itp sep.itp --nsapdb sep.pdb

This writes:
  * ``mutant.pdb``              — hybrid structure file (pass to pdb2gmx)
  * ``<FF>/<hybrid>.rtp``       — residue topology parameters (for pdb2gmx)
  * ``<FF>/<hybrid>.arn``       — atom-renaming spec (CHARMM force fields)
  * ``./<hybrid>.mtp``          — hybrid MTP file (pass to generate_hybrid_topology)

Updates ``atomtypes.atp`` and ``ffnonbonded.itp`` in the force-field directory
with any new DUM_ atom types.

**Step 2** — fill hybrid B-states in the topology (after pdb2gmx)::

    pmx generate_hybrid_topology -f topol.top -o pmxtop.top \\
        --extra_mtp A2SEP.mtp

or in Python::

    from pmx.forcefield import Topology
    from pmx.alchemy import gen_hybrid_top

    top = Topology("topol.top", ff="charmm36m-mut", version="new")
    pmxtop, itps = gen_hybrid_top(top, extra_mtp_files=["A2SEP.mtp"])

Requirements
------------
RDKit must be installed::

    conda install -c conda-forge rdkit
"""

import logging
import sys
import argparse

from pmx.model import Model
from pmx.nsaa import mutate_nsaa
from pmx.utils import get_ff_path, ff_selection
from pmx.scripts.cli import check_unknown_cmd


# =================
# Argument parsing
# =================

def parse_options():
    parser = argparse.ArgumentParser(
        description=(
            "Build a hybrid structure file for a non-standard amino acid (NSAA / PTM) "
            "free energy calculation using the pmx dual-topology approach.\n\n"
            "RDKit is required (conda install -c conda-forge rdkit)."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        '-f',
        metavar='infile',
        dest='infile',
        type=str,
        default='protein.pdb',
        help='Input structure file in PDB or GRO format. Default: protein.pdb',
    )
    parser.add_argument(
        '-o',
        metavar='outfile',
        dest='outfile',
        type=str,
        default='mutant.pdb',
        help='Output hybrid structure file in PDB format. Default: mutant.pdb',
    )
    parser.add_argument(
        '--resid',
        metavar='RESID',
        dest='resid',
        type=int,
        required=True,
        help='Residue ID to mutate (1-based after renumbering).',
    )
    parser.add_argument(
        '--resname',
        metavar='RESNAME',
        dest='resname',
        type=str,
        required=True,
        help=(
            "Three-letter name of the target NSAA (e.g. SEP, TPO, PTR).\n"
            "Must match the residue name in the supplied --itp and --nsapdb files."
        ),
    )
    parser.add_argument(
        '-ff',
        metavar='ff',
        dest='ff',
        type=str.lower,
        default=None,
        help=(
            "Force field to use (e.g. charmm36m-mut).\n"
            "If omitted, an interactive selection menu is shown."
        ),
    )
    parser.add_argument(
        '--itp',
        metavar='ITP',
        dest='itp',
        type=str,
        default=None,
        help=(
            "GROMACS ITP fragment file for the NSAA, parameterized with a\n"
            "CH3 cap at the CA position.  Required when --nsapdb is provided\n"
            "and the residue is not in the pmx built-in library."
        ),
    )
    parser.add_argument(
        '--nsapdb',
        metavar='NSAPDB',
        dest='nsapdb',
        type=str,
        default=None,
        help=(
            "PDB or GRO file containing the NSAA structure.  If omitted, the\n"
            "NSAA is taken from the pmx built-in library (standard residues only)."
        ),
    )
    parser.add_argument(
        '--chain',
        metavar='CHAIN',
        dest='chain',
        type=str,
        default=None,
        help=(
            "Chain ID of the residue to mutate.  Required only when the\n"
            "structure contains multiple chains and residue IDs are not unique."
        ),
    )
    parser.add_argument(
        '--keep_resid',
        dest='renumber',
        action='store_false',
        default=True,
        help=(
            "Do NOT renumber residues before mutation.  By default, residues\n"
            "are renumbered starting from 1 to ensure unique IDs."
        ),
    )
    parser.add_argument(
        '-v', '--verbose',
        dest='verbose',
        action='store_true',
        default=False,
        help='Print extra information during mutation.',
    )

    args, unknown = parser.parse_known_args()
    check_unknown_cmd(unknown)
    return args


# ============
# Main routine
# ============

def main(args):
    # ── Force-field selection ─────────────────────────────────────────────────
    ff = args.ff
    if ff is None:
        ff = ff_selection()

    # Verify the force field exists; get_ff_path raises if it does not
    get_ff_path(ff)

    # ── Logging ───────────────────────────────────────────────────────────────
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(levelname)s %(name)s: %(message)s',
        stream=sys.stdout,
    )

    # ── Load input structure ──────────────────────────────────────────────────
    print('log_> Loading structure: %s' % args.infile)
    m = Model(args.infile,
              renumber_residues=args.renumber,
              bPDBTER=True,
              rename_atoms=True,
              scale_coords='A')

    if args.verbose:
        print('log_> Structure loaded: %d residues, %d atoms'
              % (len(m.residues), len(m.atoms)))

    # ── Run the NSAA mutation ─────────────────────────────────────────────────
    print('log_> Mutating residue %d to %s ...' % (args.resid, args.resname))
    print(args.chain)
    print(args.chain)
    print(args.chain)
    print(args.chain)

    mutate_nsaa(
        m=m,
        mut_resid=args.resid,
        mut_resname=args.resname,
        ff=ff,
        itp=args.itp,
        nonstandardPDB=args.nsapdb,
        mut_chain=args.chain,
        inplace=True,
        verbose=args.verbose,
    )

    # ── Write output ──────────────────────────────────────────────────────────
    m.write(args.outfile)
    print('')
    print('log_> Hybrid structure written to: %s' % args.outfile)
    print('log_> Force-field files updated in: %s' % get_ff_path(ff))
    print('log_> MTP file written to working directory (pass to '
          'generate_hybrid_topology via --extra_mtp).')
    print('')


def entry_point():
    args = parse_options()
    main(args)


if __name__ == '__main__':
    entry_point()
