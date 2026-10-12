#!/usr/bin/env python

"""This module contains functions that can be used to setup the hybrid
structures and topologies needed for alchemical free energy calculations.
"""

import os
import sys
from copy import deepcopy
import numpy as np
from . import library
from .atom import Atom
from .model import Model
from .utils import mtpError, UnknownResidueError, MissingTopolParamError
from .geometry import Rotation, nuc_super, bb_super
from .mutdb import read_mtp_entry
from .forcefield import Topology, _check_case, _atoms_morphe
from .utils import get_mtp_file

__all__ = ['mutate', 'gen_hybrid_top', 'write_split_top', 'AbsRestraints']

_NOS_CYS_HYBRIDS = {'CNOS', 'U2UO'}   # NOS-Cys (carries SG2 + bridging OD2)
_NOS_LYS_HYBRIDS = {'KNOS', 'K2KN'}   # NOS-Lys (carries NZ2)

# ==============
# Main Functions
# ==============
# Cross-link and Cys/Tyr special codes that stay on the pre-built library path.
# NOTE: CD/DC (disulfide) and CN/KN (NOS) are cross-links -- they need the
# decouple+bond composition (disulfide_edit/nos_edit), not the single-residue
# `mutate` that covalent.mutate performs.  That generative cross-link path exists
# in pmx.covalent but is not wired to the CLI yet, so these still use the shipped
# split-atom hybrids (C2CD/C2CN/K2KN).  This is temporary: once the cross-link
# recipes are wired to `pmx mutate`, these become generative too and the entries
# retire -- the split here is "what's wired", not a classic-vs-generative divide.
_CLASSIC_SPECIAL = {'CM', 'CD', 'DC', 'YM', 'CN', 'KN'}


def _is_classic_mutation(mut_resname):
    """True if *mut_resname* is a standard AA / protonation / cross-link code the
    pre-built library path handles; False for PTMs and arbitrary NSAA targets,
    which the generative engine (covalent.mutate) builds on the fly."""
    aa = mut_resname.upper()
    return (aa in _CLASSIC_SPECIAL or len(aa) == 1
            or aa in library._ext_one_letter)


def mutate(m, mut_resid, mut_resname, ff, mut_chain=None,
           refB=None, inplace=False, verbose=False, itp=None, structure=None,
           generative=False, reverse=False, hybrid_name=None):
    """Creates an hybrid structure file.

    Parameters
    ----------
    m : Model
        the model to be mutated. See :py:class:`pmx.model.Model`.
    mut_resid : int
        the ID of the residue to be mutated.
    mut_resname : str
        the name of the target residue.
    ff : str
        the forcefield to use.
    mut_chain : str, optional
        the chain ID for the residue you want to mutate. This is needed if you
        have multiple chains and are providing a Model that did not have its
        residues renumbered.
    refB : str, optional
        reference structure file of the B state in PDB or GRO format. If it is
        provided, the dummy atoms will be built based on the position of the
        atoms in this structure. This option is available only for Protein
        mutations.
    inplace : bool, optional
        whether to modify the input Model in place. Default is False.
    verbose : bool, optional
        whether to print info. Default is False.

    Returns
    -------
    m2 : Model
        model with mutated residues. If argument 'inplace' is True, nothing is
        returned.
    """

    if inplace is True:
        m2 = m
    elif inplace is False:
        m2 = deepcopy(m)

    # convert Model to Angstrom units: this code works on A units
    if m2.unity == 'nm':
        m2.nm2a()
        if verbose:
            print('Model units have been converted from nm to A.')
    # get the residue based on the index
    # fetch_residue also checks the selection is valid and unique
    residue = m2.fetch_residue(idx=mut_resid, chain=mut_chain)

    # get the correct mtp file
    mtp_file = get_mtp_file(residue, ff)

    # Molecule.moltype comes from the hardcoded name sets in pmx.library, so a
    # residue the generative engine built is 'unknown' there.  Ask the force field
    # instead -- exactly as get_mtp_file does -- or the dispatch below matches no
    # branch at all.  That used to fall straight through and return None having
    # changed nothing: ABA53->ATB "succeeded" and wrote back a pure ABA residue,
    # so twelve GB1 all-to-all edges were silent zero-dG no-ops.
    moltype = residue.moltype
    if moltype == 'unknown':
        from .utils import residuetype_from_ff, get_ff_path
        moltype = residuetype_from_ff(residue.resname, get_ff_path(ff=ff)) or moltype

    # Terminal residue deletion -- a residue morph through the generative core:
    # the neighbour morphs internal->terminus (an ordinary build_hybrid_residue
    # graph map, emitted into nsaa.{rtp,mtp}), the target is dummified at
    # gen_hybrid_top.  No pre-built XdeC/XdeN library, no bespoke structure prep.
    if mut_resname.upper() == 'DEL':
        if moltype != 'protein':
            raise RuntimeError("Terminal deletion is only supported for "
                               "protein residues")
        if reverse:
            # delete_terminal has no insertion direction: the target is
            # dummified in state B unconditionally.  Accepting the flag and
            # dropping it would silently leave both legs of a double-system
            # single-box deleting, so the box charge would not balance.
            raise ValueError(
                "--reverse is not supported for a terminal deletion: there is "
                "no insertion direction, so the reference leg of a "
                "double-system-single-box cannot be built this way. Only a "
                "deletion that removes a CHARGED side chain needs a balancing "
                "reverse leg; for one, run a single box with charge-scaled "
                "(ECC) ions instead.")
        from .covalent import delete_terminal
        delete_terminal(m2, del_resid=mut_resid, ff=ff, del_chain=mut_chain,
                        inplace=True, verbose=verbose)
    # Mutation if Protein
    elif moltype == 'protein':
        # PTM aliases (P1/P2, me1/2/3), a direct PTM/NSAA target name, or an
        # --itp fragment go through the generative engine (covalent.mutate),
        # which builds the hybrid on the fly and accumulates it into nsaa.{rtp,
        # mtp}.  Standard mutations, protonation changes and cross-links (CD/CN/
        # KN) keep the pre-built library path.
        from .covalent import ptm
        # A residue at a terminus also goes through the generative engine: it
        # patches both states with the terminus (from_rtp(term=...)) and maps by
        # graph, whereas the pre-built internal X2Y hybrids assume the internal
        # backbone (HN, O) and fail on the terminal H1/H2/H3 / OT1/OT2.
        chain = residue.chain
        is_terminal = (residue == chain.cterminus()
                       or residue == chain.nterminus())
        # Detection funnel: a PTM alias, an --itp fragment, a target not shipped
        # as a pre-built 20x20 hybrid, a terminal residue, or an explicit
        # --generative override all go through the generative engine; otherwise
        # use the shipped hybrid.
        target = ptm.target_for(residue.resname, mut_resname)
        # The shipped 20x20 library is indexed by the SOURCE's one-letter code, so
        # it cannot serve a non-canonical source however ordinary the target is:
        # ABA->LEU used to reach _convert_aa_name('ABA') and die on KeyError.
        source_is_classic = library._canonical_code(residue.resname) is not None
        if target is None and (itp is not None or generative or reverse
                               or is_terminal
                               or not source_is_classic
                               or not _is_classic_mutation(mut_resname)):
            target = mut_resname
        if target is not None:
            from .covalent import mutate as _generative_mutate
            _generative_mutate(m2, mut_resid=mut_resid, target=target, ff=ff,
                               structure=structure, itp=itp, refB=refB,
                               mut_chain=mut_chain, inplace=True, verbose=verbose,
                               reverse=reverse, hybrid_name=hybrid_name)
        else:
            new_aa_name = _convert_aa_name(mut_resname)
            apply_aa_mutation(m=m2, residue=residue, new_aa_name=new_aa_name,
                              mtp_file=mtp_file, refB=refB, verbose=verbose)
    # Mutation if DNA or RNA
    elif moltype in ['dna', 'rna']:
        new_nuc_name = mut_resname.upper()
        apply_nuc_mutation(m=m2, residue=residue, new_nuc_name=new_nuc_name,
                           mtp_file=mtp_file, verbose=verbose)
    else:
        # Never fall through silently.  Doing nothing and returning the model
        # unchanged is indistinguishable from success to every caller.
        raise ValueError(
            "Cannot mutate residue %s%s: pmx classifies it as moltype %r and the "
            "force field does not identify it as protein, DNA or RNA either."
            % (residue.resname, residue.id, residue.moltype))

    if inplace is False:
        return m2


def apply_aa_mutation(m, residue, new_aa_name, mtp_file, refB=None,
                      verbose=False):

    if residue.resname == 'ILE':
        _rename_ile(residue)
    olkey = _convert_aa_name(residue.resname)

    # olkey should contain the correct one letter name of the WT residue
    # however, due to the different namings of the residues in the FFs
    # Lys needs to be checked once again: in OPLS Lys is non-protonated,
    # while in the other FFs it is protonated
    if ('opls' in mtp_file) and ('LYS' in residue.resname):
        olkey = _check_OPLS_LYS(residue)

    hybrid_residue_name = olkey+'2'+new_aa_name
    if verbose is True:
        print('log_> Residue to mutate: %s | %s | %s ' % (str(residue.id), residue.resname, residue.chain_id))
        print('log_> Mutation to apply: %s->%s' % (olkey, new_aa_name))
        print('log_> Hybrid residue name: %s' % hybrid_residue_name)
    hybrid_res, bonds, imps, diheds, rotdic = _get_hybrid_residue(residue_name=hybrid_residue_name,
                                                                  mtp_file=mtp_file,
                                                                  version='new',
                                                                  verbose=False)
    bb_super(residue, hybrid_res)

    # rename residue atoms
    hash1 = _rename_to_match_library(residue)
    hash2 = _rename_to_match_library(hybrid_res)
    _set_conformation(residue, hybrid_res, rotdic)
    if refB is not None:
        if verbose is True:
            print("log_> Set Bstate geometry according to the provided structure")
        mB = Model(refB, bPDBTER=True, rename_atoms=True, scale_coords='A')
        residueB = mB.residues[residue.id-1]
        bb_super(residue, residueB)
        for atom in hybrid_res.atoms:
            if (atom.name[0] == 'D') or atom.name.startswith('HV'):
                for atomB in residueB.atoms:
                    if atomB.name == hybrid_res.morphes[atom.name]['n1']:
                        atom.x = atomB.x
    _rename_back(residue, hash1)
    _rename_back(hybrid_res, hash2)

    # do the actual replacement
    m.replace_residue(residue=residue, new=hybrid_res, bKeepResNum=True)

    if verbose is True:
        print('log_> Inserted hybrid residue %s at position %s (chain %s)' %
              (hybrid_res.resname, str(hybrid_res.id), hybrid_res.chain_id))


def apply_nuc_mutation(m, residue, new_nuc_name, mtp_file, verbose=False):

    hybrid_residue_name, resname1, resname2 = get_nuc_hybrid_resname(residue, new_nuc_name)
    if verbose is True:
        print('log_> Residue to mutate: %d | %s | %s ' % (residue.id, residue.resname, residue.chain_id))
        print('log_> Mutation to apply: %s->%s' % (residue.resname[1], new_nuc_name))
        print('log_> Hybrid residue name: %s' % hybrid_residue_name)
    hybrid_res, bonds, imps, diheds, rotdic = _get_hybrid_residue(residue_name=hybrid_residue_name,
                                                                  mtp_file=mtp_file,
                                                                  version='new',
                                                                  verbose=False)

    nuc_super(residue, hybrid_res, resname1, resname2)
    for atom in hybrid_res.atoms:
        if atom.name[0] != 'D':
            atom.x = residue[atom.name].x
    m.replace_residue(residue=residue, new=hybrid_res, bKeepResNum=True)
    if verbose is True:
        print('log_> Inserted hybrid residue %s at position %d (chain %s)' %
              (hybrid_res.resname, hybrid_res.id, hybrid_res.chain_id))


def get_nuc_hybrid_resname(residue, new_nuc_name):

    if residue.moltype == 'dna':
        firstLetter = 'D'
    if residue.moltype == 'rna':
        firstLetter = 'R'

    # identify if the nucleotide is terminal
    for a in residue.atoms:
        if a.name == 'H3T':
            r1 = firstLetter+residue.resname[1]+'3'
            r2 = firstLetter+new_nuc_name+'3'
            dict_key = r1+'_'+r2
            if residue.moltype == 'rna':
                hybrid_residue_name = _rna_names[dict_key]
            else:
                hybrid_residue_name = _dna_names[dict_key]
            return(hybrid_residue_name, residue.resname[1], new_nuc_name)
        elif a.name == 'H5T':
            r1 = firstLetter+residue.resname[1]+'5'
            r2 = firstLetter+new_nuc_name+'5'
            dict_key = r1+'_'+r2
            if residue.moltype == 'rna':
                hybrid_residue_name = _rna_names[dict_key]
            else:
                hybrid_residue_name = _dna_names[dict_key]
            return(hybrid_residue_name, residue.resname[1], new_nuc_name)
    hybrid_residue_name = residue.resname+new_nuc_name
    return(hybrid_residue_name, residue.resname[1], new_nuc_name)


def gen_hybrid_top(topol, recursive=True, verbose=False, scaleDih=1.0,
                   extra_mtp_files=None, supplement_bonded_files=None,
                   inter_bonds=None, decouple_dummies=False):
    """Fills the bstate of a topology file containing pmx hybrid residues. This
    can be either a top or itp file. If the file contains other itp files via
    include statements, the function can iterate through them if the recursive
    flag is set to True.

    It returns a Topology instance with b states filled for the input topol
    instance, and optionally a list of additional Topology instances for the
    included itp files. If there are no itp files, or if recursive is set to
    False, then an empty list is returned.

    Parameters
    ----------
    topol : Topology
        input Topology instance
    recursive : bool, optional
        whether to fill b states also for all itp files included in the
        topology file. Default is True.
    verbose : bool
        whether to print out info. Default is False.
    scaleDih : float
        scaling for dihedrals with a dummy.
    extra_mtp_files : list of str, optional
        Additional MTP files to search for residue parameters, in addition to
        the force-field's built-in MTP database.  These files take precedence
        over the built-in database.  Pass the ``.mtp`` files generated by
        :func:`~pmx.nsaa.mutate_nsaa` here so that NSAA hybrid topologies are
        filled correctly::

            pmxtop, itps = gen_hybrid_top(top, extra_mtp_files=["A2SEP.mtp"])

    decouple_dummies : bool, optional
        Opt-in Boresch dummy-atom decoupling corrections (default False).  When
        True, removes factorization-breaking bonded terms at physical/dummy
        junctions that the default "naive" symmetric dummy treatment leaves in
        place (see :func:`_decouple_dummies`).  Off by default so existing
        behaviour is byte-identical.

    Returns
    -------
    top, itps : obj, list of obj
        output Topology instance with b states assigned (top), and
        list of Topology instances for the additional itp files (itps)
        included in the input topology.
    """

    # ff setup
    ff = topol.forcefield
    ffbonded_file = os.path.join(topol.ffpath, 'ffbonded.itp')

    # main function that returns the Topology with filled B states
    def process_topol(topol, ff, ffbonded_file, verbose=False, scaleDih=1.0,
                      extra_mtp_files=None, supplement_bonded_files=None,
                      decouple_dummies=False):
        pmxtop = deepcopy(topol)
        # create model with residue list
        m = Model(atoms=pmxtop.atoms, renumber_residues=False)
        # need to sequentially renumber residues,
        # because --keep_resid flag may leave residues numbered non-sequentially
        counter = 0
        for r in m.residues:
            counter += 1
            r.id = counter
        # use model residue list
        pmxtop.residues = m.residues
        # get list of hybrid residues and their params
        rlist, rdic = _get_hybrid_residues(m=m, ff=ff, version='new',
                                           verbose=verbose,
                                           extra_mtp_files=extra_mtp_files)
        # correct b-states
        pmxtop.assign_fftypes()
        # inject NSAA supplement bonded parameters before B-state lookup
        if supplement_bonded_files:
            _inject_supplement_bonded(pmxtop, supplement_bonded_files)
        if verbose is True:
            for r in rlist:
                print('log_> Hybrid Residue -> %d | %s ' % (r.id, r.resname))

        # Dummify terminal deletion target residues BEFORE bonded entry lookup.
        # When a deletion neighbor (XdeN/XdeC) is present, the flanking deletion
        # target residue (e.g. Met1 for an N-terminal deletion) must have its
        # B-state atom types set to DUM_ first.  Otherwise _find_bonded_entries
        # sees atomtypeB=None for the deleted residue's atoms and cannot resolve
        # the inter-residue bond (e.g. Met1.C -- EdeN.N).
        for res in rlist:
            if res.resname.endswith('deC'):
                chain_res = res.chain.residues
                idx = chain_res.index(res)
                if idx + 1 < len(chain_res):
                    del_res = chain_res[idx + 1]
                    if verbose:
                        print('log_> Dummifying deletion target: %d | %s'
                              % (del_res.id, del_res.resname))
                    for atom in del_res.atoms:
                        atom.atomtypeB = 'DUM_' + atom.atomtype
                        atom.typeB = 'DUM_' + atom.type
                        atom.qB = 0.0
                        atom.mB = atom.m
            elif res.resname.endswith('deN'):
                chain_res = res.chain.residues
                idx = chain_res.index(res)
                if idx - 1 >= 0:
                    del_res = chain_res[idx - 1]
                    if verbose:
                        print('log_> Dummifying deletion target: %d | %s'
                              % (del_res.id, del_res.resname))
                    for atom in del_res.atoms:
                        atom.atomtypeB = 'DUM_' + atom.atomtype
                        atom.typeB = 'DUM_' + atom.type
                        atom.qB = 0.0
                        atom.mB = atom.m

        _find_bonded_entries(pmxtop, verbose=verbose)
        _find_angle_entries(pmxtop, verbose=verbose)
        dih_predef_default = []
        _find_predefined_dihedrals(pmxtop, rlist, rdic, ffbonded_file,
                                   dih_predef_default, ff, scaleDih=scaleDih)
        _find_dihedral_entries(pmxtop, rlist, rdic, dih_predef_default,
                               verbose=verbose, scaleDih=scaleDih)

        _add_extra_DNA_RNA_impropers(pmxtop, rlist, 1, [180, 40, 2], [180, 40, 2])

        if verbose is True:
            print('log_> Total charge of state A = %.f' % pmxtop.get_qA())
            print('log_> Total charge of state B = %.f' % pmxtop.get_qB())

        # fix B-states for terminus hybrid residues
        for res in rlist:
            nterm = cterm = False
            if res == res.chain.nterminus():
                nterm = True
            elif res == res.chain.cterminus():
                cterm = True
            else:
                continue
            if cterm:
                for atom in res.atoms:
                    if atom.name in ("C", "OT1", "OT2", "OC1", "OC2"):
                        atom.typeB = atom.type
                        atom.atomtypeB = atom.atomtype
                        atom.mB = atom.m
                        atom.qB = atom.q
            if nterm:
                for atom in res.atoms:
                    if atom.name in ("N", "H1", "H2", "H3"):
                        atom.typeB = atom.type
                        atom.atomtypeB = atom.atomtype
                        atom.mB = atom.m
                        atom.qB = atom.q
                    if atom.name in ("CA", "HA"):
                        atom.typeB = atom.type
                        atom.atomtypeB = atom.atomtype

        # if prolines are involved, break one bond (CD-CG)
        # and angles X-CD-CG, CD-CG-X
        # and dihedrals with CD and CG
        _proline_decouplings(pmxtop, rlist, rdic)
        # also decouple all dihedrals with [CD and N] and [CB and Calpha] for proline
        _proline_dihedral_decouplings(pmxtop, rlist, rdic)

        # Inter-residue cross-link bonds (disulfide, NOS, SONOS, ...) are formed
        # from the user's gentop --bond declarations through the one graph
        # consumer below -- no resname auto-detection.  Terminal deletion is a
        # structural decoupling (depends on the runtime chain, not a named-atom
        # bond), so it keeps its own detector, which already delegates to the
        # unified apply_bond_morph(mode='decouple').
        _terminal_deletion_decouplings(pmxtop, rlist, verbose=verbose)
        _apply_inter_bonds(pmxtop, rlist, inter_bonds, verbose=verbose)

        # Opt-in Boresch dummy-atom decoupling (off by default).  Removes
        # factorization-breaking bonded terms at physical/dummy junctions that
        # the default "naive" symmetric treatment leaves in place.
        if decouple_dummies:
            _decouple_dummies(pmxtop, verbose=verbose)
        return pmxtop

    # -------------
    # Main topology
    # -------------
    if verbose is True:
        print('\nlog_> Reading input %s file "%s"'
              % (topol.filename.split('.')[-1], topol.filename))
    pmx_top = process_topol(topol=topol, ff=ff, ffbonded_file=ffbonded_file,
                            verbose=verbose, scaleDih=scaleDih,
                            extra_mtp_files=extra_mtp_files,
                            supplement_bonded_files=supplement_bonded_files,
                            decouple_dummies=decouple_dummies)

    # ------------------------------------
    # itps too if asked for and if present
    # ------------------------------------
    itps = topol.include_itps
    pmx_itps = []
    if recursive is True and len(itps) > 0:
        for itp, where in itps:
            if verbose is True:
                print('\nlog_> Reading input itp file "%s""' % (itp))
            itp_path = os.path.dirname(os.path.relpath(topol.filename))
            itp_filename = os.path.join(itp_path, itp)
            topol2 = Topology(itp_filename, ff=ff, version='new')
            # fill b states
            pmx_itp = process_topol(topol=topol2, ff=ff,
                                    ffbonded_file=ffbonded_file,
                                    verbose=verbose,
                                    extra_mtp_files=extra_mtp_files,
                                    supplement_bonded_files=supplement_bonded_files,
                                    decouple_dummies=decouple_dummies)
            pmx_itps.append(pmx_itp)

    # A --bond declaration is applied by whichever molecule owns both endpoints
    # (see _apply_inter_bonds); flag any that no molecule/chain resolved (a typo
    # or a wrong resid/chain), which per-molecule skipping would otherwise hide.
    unresolved = [eb for eb in (inter_bonds or []) if not eb.get('_applied')]
    if unresolved:
        raise RuntimeError(
            "gentop --bond: no molecule/chain contains both atoms of: "
            + "; ".join("%s:%s-%s:%s" % (eb['resid_a'], eb['atom_a'],
                                         eb['resid_b'], eb['atom_b'])
                        for eb in unresolved))

    return pmx_top, pmx_itps


def write_split_top(pmxtop, outfile='pmxtop.top', scale_mass=False,
                    verbose=False):
    """
    Write three topology files to be used for three separate free energy
    calculations: charges off, vdw on, changes on. This can be useful when one
    wants to avoid using a soft-core for the electrostatic interactions.

    Parameters
    ----------
    pmxtop : Topology
        Topology object with B states defined.
    outfile : str
        name of output topology file.
    scale_mass : bool
        whether to scale the masses of morphing atoms
    verbose : bool
        whether to print information or not

    Returns
    -------
    None
    """

    root, ext = os.path.splitext(outfile)
    out_file_qoff = root + '_qoff' + ext
    out_file_vdw = root + '_vdw' + ext
    out_file_qon = root + '_qon' + ext

    # get charge of states A/B for hybrid residues
    qA = pmxtop.get_hybrid_qA()
    qB = pmxtop.get_hybrid_qB()

    if verbose is True:
        print('------------------------------------------------------')
        print('log_> Creating splitted topologies............')

        print('log_> Making "qoff" topology : "%s"' % out_file_qoff)
    contQ = deepcopy(qA)
    pmxtop.write(out_file_qoff, stateQ='AB', stateTypes='AA', dummy_qB='off',
                 scale_mass=scale_mass, target_qB=qA, stateBonded='AA',
                 full_morphe=False)
    if verbose is True:
        print('log_> Charge of state A: %g' % pmxtop.qA)
        print('log_> Charge of state B: %g' % pmxtop.qB)
    if verbose is True:
        print('------------------------------------------------------')
        print('log_> Making "vdw" topology : "%s"' % out_file_vdw)
    contQ = deepcopy(qA)
    pmxtop.write(out_file_vdw, stateQ='BB', stateTypes='AB', dummy_qA='off',
                 dummy_qB='off', scale_mass=scale_mass,
                 target_qB=contQ, stateBonded='AB', full_morphe=False)
    if verbose is True:
        print('log_> Charge of state A: %g' % pmxtop.qA)
        print('log_> Charge of state B: %g' % pmxtop.qB)
        print('------------------------------------------------------')

        print('log_> Making "qon" topology : "%s"' % out_file_qon)
    pmxtop.write(out_file_qon, stateQ='BB', stateTypes='BB', dummy_qA='off',
                 dummy_qB='on', scale_mass=scale_mass, target_qB=qB,
                 stateBonded='BB', full_morphe=False)
    if verbose is True:
        print('log_> Charge of state A: %g' % pmxtop.qA)
        print('log_> Charge of state B: %g' % pmxtop.qB)
        print('------------------------------------------------------')


def decouple_mol(top):
    """Setup topology for decoupling a small molecule. The nonbonded
    interactions are on in state A and are off in state B.

    Params
    ------
    top : Topology
        topology object that will be modified

    Returns
    -------
    None
    """

    for atom in top.atoms:
        atom.qB = 0.0
        atom.mB = atom.m
        atom.typeB = 'dum'
        atom.atomtypeB = 'dum'


    # create the dummy atomtype
    dummy = {}
    dummy['name'] = 'dum'
    dummy['bond_type'] = 'dum'
    dummy['mass'] = 0.0
    dummy['charge'] = 0.0
    dummy['ptype'] = 'A'
    dummy['sigma'] = 0.0
    dummy['epsilon'] = 0.0

    # add the dummy to the topology atomtypes
    top.atomtypes.append(dummy)


def couple_mol(top):
    """Setup topology for coupling a small molecule. The nonbonded
    interactions are off in state A and are on in state B.

    Params
    ------
    top : Topology
        topology object that will be modified

    Returns
    -------
    None
    """

    for atom in top.atoms:
        atom.qB = atom.q
        atom.q = 0.0
        atom.mB = atom.m
        atom.typeB = atom.type
        atom.atomtypeB = atom.type
        atom.type = 'dum'
        atom.atomtype = 'dum'

    # create the dummy atomtype
    dummy = {}
    dummy['name'] = 'dum'
    dummy['bond_type'] = 'dum'
    dummy['mass'] = 0.0
    dummy['charge'] = 0.0
    dummy['ptype'] = 'A'
    dummy['sigma'] = 0.0
    dummy['epsilon'] = 0.0

    # add the dummy to the topology atomtypes
    top.atomtypes.append(dummy)


# ===============
# HelperFunctions
# ===============
def _scale_dih( d, scale ):
    # scale dihedral
    dih = deepcopy(d)
    if dih[0]==9 or dih[0]==4 or dih[0]==1 or dih[0]==2:
        dih[2] *= scale
    elif dih[0]==3:
        dih[1] *= scale
        dih[2] *= scale
        dih[3] *= scale
        dih[4] *= scale
        dih[5] *= scale
        dih[6] *= scale
    return(dih)

def _check_OPLS_LYS(res):
    if res.has_atom('HZ3'):
        return('K')
    else:
        return('O')


def _convert_aa_name(aa):
    # special multi-character pmx codes that are not in the standard AA library
    # CM : deprotonated cysteine (CYM) target for C2CM hybrid
    # CD : disulfide cysteine target for C2CD hybrid (alchemical SS formation)
    # DC : reverse disulfide target for D2DC hybrid (alchemical SS breaking)
    # P1 : monoanionic phospho-amino acid target for S2P1 (SER), T2P1 (THR), or Y2P1 (TYR)
    # P2 : dianionic phospho-amino acid target for S2P2 (SER), T2P2 (THR), or Y2P2 (TYR)
    # YM: deprotonated tyrosine
    _special = {'CM', 'CD', 'DC', 'P1', 'P2', 'YM', 'CN', 'KN'}
    if aa.upper() in _special:
        return aa.upper()
    elif len(aa) == 1:
        return aa.upper()
    elif len(aa) in [3, 4]:
        return library._ext_one_letter[aa.upper()]
    else:
        raise UnknownResidueError(aa)


def _rename_ile(residue):
    dic = {'CD':  'CD1',
           'HD1': 'HD11',
           'HD2': 'HD12',
           'HD3': 'HD13'
           }
    for key, value in dic.items():
        try:
            atom = residue[key]
            atom.name = value
        except:
            pass


def _rename_to_match_library(res):
    name_hash = {}
    atoms = res.atoms
    for atom in atoms:
        foo = atom.name
        # for serine
        if (atom.resname == 'SER') and (atom.name == 'HG1'):
            atom.name = 'HG'
        if ('S2' in atom.resname) and (atom.name == 'HG1'):
            atom.name = 'HG'
        # for cysteine
        if (atom.resname == 'CYS') and (atom.name == 'HG1'):
            atom.name = 'HG'
        if ('C2' in atom.resname) and (atom.name == 'HG1'):
            atom.name = 'HG'

        name_hash[atom.name] = foo
    return name_hash


def _rename_back(res, name_hash):
    for atom in res.atoms:
        atom.name = name_hash[atom.name]


# Mapping from split-atom hybrid names to source-residue atom names.
# Used in _set_conformation for C2CD/D2DC (CHARMM disulfide hybrids where
# CB and SG are split into CB1/CB2 and SG1/SG2 branches).
_SPLIT_ATOM_NAME_MAP = {
    'CB1': 'CB',  'HB11': 'HB1', 'HB12': 'HB2',
    'CB2': 'CB',  'HB21': 'HB1', 'HB22': 'HB2',
    'CE1': 'CE',  'HE11': 'HE1', 'HE12': 'HE2',
    'NZ1': 'NZ',  'HZ11': 'HZ1', 'HZ12': 'HZ2', 'HZ13': 'HZ3',
    'SG1': 'SG',  'SG2':  'SG',
}


def _set_conformation(old_res, new_res, rotdic):
    #print(new_res.name)
    old_res.get_real_resname()
    dihedrals = library._aa_dihedrals[old_res.real_resname]
    for key, lst in rotdic.items():
        new = new_res.fetchm(lst)
        rotdic[key] = new
    chis = []
    for key in rotdic.keys():
        at1, at2 = key.split('-')
        for d in dihedrals:
            if d[1] == at1 and d[2] == at2 \
                   and d[-1] != -1:
                chis.append(d)
    for d in chis:
        atoms = old_res.fetchm(d[:4])
        phi = atoms[0].dihedral(atoms[1], atoms[2], atoms[3])
        atoms2 = new_res.fetchm(d[:4])
        phi2 = atoms2[0].dihedral(atoms2[1], atoms2[2], atoms2[3])
        diff = phi-phi2
        a1, a2 = new_res.fetchm(d[1:3])
        key = a1.name+'-'+a2.name
        atoms = rotdic[key]
        rot = Rotation(a1.x, a2.x)
        for atom in atoms:
            atom.x = rot.apply(atom.x, diff)
    for atom in new_res.atoms:
        # Skip atoms that are B-state-only dummies in state A
        if atom.atomtype.startswith('DUM_'):
            continue
        # Skip atoms already placed by the rotation algorithm (D/HV prefix)
        if (atom.name[0] == 'D') or atom.name.startswith('HV'):
            continue
        # Copy position from old_res; handle split-atom name mapping.
        # C2CD/D2DC (disulfide) and the NOS hybrids all split CB/SG (and CE/NZ)
        # into duplicated branches whose coordinates come from the single WT atom.
        if new_res.resname in ["C2CD", "D2DC",
                               "O2CN", "CN2O", "C2CN", "CN2C",
                               "O2KN", "KN2O", "K2KN", "KN2K"]:
            src_name = _SPLIT_ATOM_NAME_MAP.get(atom.name, atom.name)
        else:
            src_name = atom.name
        try:
            atom.x = old_res[src_name].x
        except (KeyError, AttributeError, IndexError):
            pass  # keep MTP/bb_super position if source atom not found


def _proline_ring_state(r, topol):
    """Return which state ('A'/'B') has a proline ring at residue *r*, else None.

    Proline is identified by its ring-closure bond N-Cdelta (the side chain closes
    a 5-membered ring back onto the backbone N) -- unique to proline, so detection
    is independent of the hybrid *name*.  This fires for the classic hybrid ('X2P'),
    the generative hybrid ('X2PR'), and terminal proline alike, and -- unlike the
    old ``'P' in resname and len==3`` gate -- never mis-fires on phospho hybrids
    (``S2P1`` etc.).  'A' = real proline ring in state A; 'B' = appearing (dummy)
    proline ring in state B.
    """
    byname = {}
    for a in r.atoms:
        byname.setdefault(a.name, a)
    n = byname.get('N')
    if n is None:
        return None

    def ring_cd(cd_name):
        cd = byname.get(cd_name)
        if cd is None:
            return None
        for b in topol.bonds:
            if {b[0].id, b[1].id} == {n.id, cd.id}:
                return cd
        return None

    def _morphs(atom):
        # The ring decoupling (break CD-CG + ring angles/dihedrals in one state) is only
        # valid when the ring APPEARS or DISAPPEARS -- i.e. its atoms are dummy in exactly
        # one state.  A Pro->Pro morph (N-terminal-proline terminal patch, XdeN deletion
        # neighbour) keeps the ring REAL in BOTH states; decoupling it there spuriously
        # opens an intact ring and destabilises the switch.  So require a real<->dummy morph.
        ta = (getattr(atom, 'atomtype', '') or '')
        tb = (getattr(atom, 'atomtypeB', '') or ta)
        return ta.upper().startswith('DUM') != tb.upper().startswith('DUM')

    cd = ring_cd('CD')       # real proline ring present in state A
    if cd is not None:
        return 'A' if _morphs(cd) else None
    dcd = ring_cd('DCD')     # appearing (dummy) proline ring in state B
    if dcd is not None:
        return 'B' if _morphs(dcd) else None
    return None


def _proline_dihedral_decouplings(topol, rlist, rdic):
    for r in rlist:
        prolineState = _proline_ring_state(r, topol)
        if prolineState is not None:
            # find the atoms for which the bond needs to be broken
            prolineCDatom = ''
            prolineCBatom = ''
            prolineNatom = ''
            prolineCAatom = ''
            for a in r.atoms:
                if (prolineState == 'A') and (a.name == 'CD'):
                    prolineCDatom = a
                elif (prolineState == 'A') and (a.name == 'N'):
                    prolineNatom = a
                elif (prolineState == 'A') and (a.name == 'CB'):
                    prolineCBatom = a
                elif (prolineState == 'A') and (a.name == 'CA'):
                    prolineCAatom = a
                elif (prolineState == 'B') and (a.name == 'DCD'):
                    prolineCDatom = a
                elif (prolineState == 'B') and (a.name == 'N'):
                    prolineNatom = a
                elif (prolineState == 'B') and (a.name == 'DCB'):
                    prolineCBatom = a
                elif (prolineState == 'B') and (a.name == 'CA'):
                    prolineCAatom = a

            # decouple [CD and N] dihedrals
            for dih in topol.dihedrals:
                a1 = dih[0]
                a2 = dih[1]
                a3 = dih[2]
                a4 = dih[3]
                func = dih[4]
                if ((a1.id == prolineCDatom.id and a2.id == prolineNatom.id)) or \
                   ((a2.id == prolineCDatom.id and a1.id == prolineNatom.id)) or \
                   ((a1.id == prolineCDatom.id and a3.id == prolineNatom.id)) or \
                   ((a3.id == prolineCDatom.id and a1.id == prolineNatom.id)) or \
                   ((a1.id == prolineCDatom.id and a4.id == prolineNatom.id)) or \
                   ((a4.id == prolineCDatom.id and a1.id == prolineNatom.id)) or \
                   ((a2.id == prolineCDatom.id and a3.id == prolineNatom.id)) or \
                   ((a3.id == prolineCDatom.id and a2.id == prolineNatom.id)) or \
                   ((a2.id == prolineCDatom.id and a4.id == prolineNatom.id)) or \
                   ((a4.id == prolineCDatom.id and a2.id == prolineNatom.id)) or \
                   ((a3.id == prolineCDatom.id and a4.id == prolineNatom.id)) or \
                   ((a4.id == prolineCDatom.id and a3.id == prolineNatom.id)):
                    if (dih[5] == 'NULL') or (dih[6] == 'NULL'):
                        continue
                    if prolineState == 'A':
                        if func == 3:
                            dih[6] = [func, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        elif func == 2:
                            dih[6] = [func, dih[6][1], 0.0]
                        else:
                            dih[6] = [func, dih[6][1], 0.0, dih[6][-1]]
                    if prolineState == 'B':
                        if func == 3:
                            dih[5] = [func, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        elif func == 2:
                            dih[5] = [func, dih[6][1], 0.0]
                        else:
                            dih[5] = [func, dih[6][1], 0.0, dih[6][-1]]

            # decouple [CB and CA] dihedrals
            for dih in topol.dihedrals:
                a1 = dih[0]
                a2 = dih[1]
                a3 = dih[2]
                a4 = dih[3]
                func = dih[4]
                if ((a1.id == prolineCBatom.id and a2.id == prolineCAatom.id)) or \
                   ((a2.id == prolineCBatom.id and a1.id == prolineCAatom.id)) or \
                   ((a1.id == prolineCBatom.id and a3.id == prolineCAatom.id)) or \
                   ((a3.id == prolineCBatom.id and a1.id == prolineCAatom.id)) or \
                   ((a1.id == prolineCBatom.id and a4.id == prolineCAatom.id)) or \
                   ((a4.id == prolineCBatom.id and a1.id == prolineCAatom.id)) or \
                   ((a2.id == prolineCBatom.id and a3.id == prolineCAatom.id)) or \
                   ((a3.id == prolineCBatom.id and a2.id == prolineCAatom.id)) or \
                   ((a2.id == prolineCBatom.id and a4.id == prolineCAatom.id)) or \
                   ((a4.id == prolineCBatom.id and a2.id == prolineCAatom.id)) or \
                   ((a3.id == prolineCBatom.id and a4.id == prolineCAatom.id)) or \
                   ((a4.id == prolineCBatom.id and a3.id == prolineCAatom.id)):
                    if (dih[5] == 'NULL') or (dih[6] == 'NULL'):
                        continue
                    if prolineState == 'A':
                        if func == 3:
                            dih[6] = [func, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        elif func == 2:
                            dih[6] = [func, dih[6][1], 0.0]
                        else:
                            dih[6] = [func, dih[6][1], 0.0, dih[6][-1]]
                    if prolineState == 'B':
                        if func == 3:
                            dih[5] = [func, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        elif func == 2:
                            dih[5] = [func, dih[6][1], 0.0]
                        else:
                            dih[5] = [func, dih[6][1], 0.0, dih[6][-1]]


def _proline_decouplings(topol, rlist, rdic):
    for r in rlist:
        prolineState = _proline_ring_state(r, topol)
        if prolineState is not None:
            # find the atoms for which the bond needs to be broken
            prolineCDatom = ''
            prolineCGatom = ''
            for a in r.atoms:
                if (prolineState == 'A') and (a.name == 'CD'):
                    prolineCDatom = a
                elif (prolineState == 'A') and (a.name == 'CG'):
                    prolineCGatom = a
                elif (prolineState == 'B') and (a.name == 'DCD'):
                    prolineCDatom = a
                elif (prolineState == 'B') and (a.name == 'DCG'):
                    prolineCGatom = a

            # break the bond
            for b in topol.bonds:
                a1 = b[0]
                a2 = b[1]
                if ((a1.id == prolineCDatom.id) and (a2.id == prolineCGatom.id)) or \
                   ((a2.id == prolineCDatom.id) and (a1.id == prolineCGatom.id)):
                    if prolineState == 'A':
                        b[4] = deepcopy(b[3])
                        b[4][2] = 0.0
                    elif prolineState == 'B':
                        b[3] = deepcopy(b[4])
                        b[3][2] = 0.0

            # decouple angles
            for ang in topol.angles:
                a1 = ang[0]
                a2 = ang[1]
                a3 = ang[2]
                func = ang[3]
                if ((a1.id == prolineCDatom.id and a2.id == prolineCGatom.id)) or \
                   ((a2.id == prolineCDatom.id and a1.id == prolineCGatom.id)) or \
                   ((a1.id == prolineCDatom.id and a3.id == prolineCGatom.id)) or \
                   ((a3.id == prolineCDatom.id and a1.id == prolineCGatom.id)) or \
                   ((a2.id == prolineCDatom.id and a3.id == prolineCGatom.id)) or \
                   ((a3.id == prolineCDatom.id and a2.id == prolineCGatom.id)):
                    if prolineState == 'A':
                        if func == 5:
                            ang[5] = [1, 0, 0, 0, 0]
                        else:
                            ang[5] = [1, 0, 0]
                    elif prolineState == 'B':
                        if func == 5:
                            ang[4] = [1, 0, 0, 0, 0]
                        else:
                            ang[4] = [1, 0, 0]

            # decouple dihedrals
            for dih in topol.dihedrals:
                a1 = dih[0]
                a2 = dih[1]
                a3 = dih[2]
                a4 = dih[3]
                func = dih[4]
                if ((a1.id == prolineCDatom.id and a2.id == prolineCGatom.id)) or \
                   ((a2.id == prolineCDatom.id and a1.id == prolineCGatom.id)) or \
                   ((a1.id == prolineCDatom.id and a3.id == prolineCGatom.id)) or \
                   ((a3.id == prolineCDatom.id and a1.id == prolineCGatom.id)) or \
                   ((a1.id == prolineCDatom.id and a4.id == prolineCGatom.id)) or \
                   ((a4.id == prolineCDatom.id and a1.id == prolineCGatom.id)) or \
                   ((a2.id == prolineCDatom.id and a3.id == prolineCGatom.id)) or \
                   ((a3.id == prolineCDatom.id and a2.id == prolineCGatom.id)) or \
                   ((a2.id == prolineCDatom.id and a4.id == prolineCGatom.id)) or \
                   ((a4.id == prolineCDatom.id and a2.id == prolineCGatom.id)) or \
                   ((a3.id == prolineCDatom.id and a4.id == prolineCGatom.id)) or \
                   ((a4.id == prolineCDatom.id and a3.id == prolineCGatom.id)):
                    if (dih[5] == 'NULL') or (dih[6] == 'NULL'):
                        continue

                    if prolineState == 'A':
                        if func == 3:
                            dih[6] = [func, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        elif func == 2:
                            dih[6] = [func, dih[6][1], 0.0]
                        else:
                            dih[6] = [func, dih[6][1], 0.0, dih[6][-1]]
                    if prolineState == 'B':
                        if func == 3:
                            dih[5] = [func, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        elif func == 2:
                            dih[5] = [func, dih[6][1], 0.0]
                        else:
                            dih[5] = [func, dih[6][1], 0.0, dih[6][-1]]


def _dummy_in_state(atom, state):
    """True if *atom* is a dummy (DUM_ atomtype) in end state 'A' or 'B'.

    Mirrors :func:`_check_case`: state A reads ``atomtype``, state B reads
    ``atomtypeB`` (falling back to ``atomtype`` when it is unset).
    """
    if state == 'A':
        t = atom.atomtype
    else:
        t = atom.atomtypeB if atom.atomtypeB is not None else atom.atomtype
    return t is not None and t.startswith('DUM')


def _decouple_dummies(topol, verbose=False):
    """Boresch dummy-atom decoupling corrections (opt-in: ``gentop --decouple-dummies``).

    pmx's default ("naive") dummy treatment keeps *every* bonded term of an
    appearing/disappearing atom in BOTH the physical and the dummy end state.
    Some of those terms couple a physical degree of freedom to a dummy one at
    the physical/dummy junction, which breaks the factorization of the dummy
    contribution to the partition function -- the condition that makes the dummy
    free energy a state-independent constant that cancels in a thermodynamic
    cycle (Fleck, Wieder, Boresch, *J. Chem. Theory Comput.* 2021, 17, 4403).
    This pass removes the offending coupling terms *in the state where the atom
    is a dummy*, keeping enough of the real bonded terms that no under-restraint
    ("flapping") is introduced (validated by gas-phase MD; see the sub-functions).

    Off by default; wired through :func:`gen_hybrid_top` and the
    ``gentop --decouple-dummies`` flag.  Two corrections are applied:

    * **Case #4** -- Urey-Bradley junction terms (:func:`_dc_urey_bradley`).
    * **Cases #2/#3** -- redundant terminal-junction dihedrals
      (:func:`_dc_terminal_dihedrals`).

    All modifications are per-state and touch only the state where the atom is a
    dummy; the physical end state is left byte-identical.
    """
    n_ub = _dc_urey_bradley(topol, verbose=verbose)
    n_del, n_anchor = _dc_terminal_dihedrals(topol, verbose=verbose)
    if verbose:
        print('log_> decouple-dummies: %d UB removed, %d redundant junction '
              'dihedrals deleted, %d anchor dihedral(s) kept'
              % (n_ub, n_del, n_anchor))
    return n_ub, n_del, n_anchor


def _dc_urey_bradley(topol, verbose=False):
    """Case #4 -- delete Urey-Bradley junction terms (CHARMM func-5 angles).

    A Urey-Bradley term is a harmonic term in the 1-3 distance of a bond angle,
    i.e. an effective virtual bond between the two *outer* atoms of the angle.
    Boresch: "the additional Urey-Bradley term of a bond angle needs to be
    deleted whenever the angle is made up from two physical and one dummy atom"
    (p. 4406); the paper deletes these even in its "naive" baseline.  A
    2-physical/1-dummy angle necessarily has the dummy as an *outer* atom (the
    central atom is physical), so its UB 1-3 pair straddles the junction.  We
    zero the UB sub-term (S0, kUB) in the offending state while keeping the
    harmonic angle (theta0, ktheta) -- the harmonic angle remains as the anchor,
    so no flapping is introduced (case #4 needs no re-anchoring).

    Note on aliasing: for a junction angle (physical in one state, dummy in the
    other) :func:`_find_angle_entries` stores the *same* list object in both the
    A-slot ``ang[4]`` and B-slot ``ang[5]``.  We therefore **reassign** the
    dummy-state slot with a fresh list rather than mutating it in place, so the
    physical state's parameters are left untouched.
    """
    n_ub = 0
    for ang in topol.angles:
        if len(ang) < 6:           # no B-state params filled -> nothing to do
            continue
        if ang[3] != 5:            # only CHARMM Urey-Bradley (func-5) angles
            continue
        a1, a2, a3 = ang[0], ang[1], ang[2]
        for state, slot in (('A', 4), ('B', 5)):
            params = ang[slot]
            if not isinstance(params, list) or len(params) < 5:
                continue
            d1 = _dummy_in_state(a1, state)
            d2 = _dummy_in_state(a2, state)
            d3 = _dummy_in_state(a3, state)
            # Case #4: exactly one dummy, on an OUTER position (central atom
            # physical) -> 2 physical + 1 dummy, UB 1-3 pair straddles junction.
            if (not d2) and (d1 != d3):
                if params[3] != 0.0 or params[4] != 0.0:
                    # reassign (do not mutate: A/B slots may be the same object)
                    ang[slot] = [params[0], params[1], params[2], 0.0, 0.0]
                    n_ub += 1
                    if verbose:
                        print('log_> decouple-dummies: removed UB (state %s) '
                              '%s-%s-%s' % (state, a1.name, a2.name, a3.name))
    return n_ub


def _dc_connected_components(idset, adj):
    """Connected components of the atom-ids in *idset* under adjacency *adj*."""
    seen = set()
    comps = []
    for s in idset:
        if s in seen:
            continue
        comp = set()
        stack = [s]
        while stack:
            x = stack.pop()
            if x in comp:
                continue
            comp.add(x)
            seen.add(x)
            stack.extend(y for y in adj.get(x, ()) if y in idset and y not in comp)
        comps.append(comp)
    return comps


def _dc_terminal_dihedrals(topol, verbose=False):
    """Cases #2/#3 -- prune redundant dihedrals at *terminal* dummy junctions.

    A *terminal junction* is a dummy group (a connected set of dummy atoms in a
    given end state) attached to the physical system through exactly ONE physical
    bridge atom, with no atom of the group bonded to a second physical atom (no
    ring back into the physical system).  For such junctions Boresch (p. 4409)
    shows one can keep *all bond angle terms* -- only the dihedrals are redundant.
    pmx's naive treatment keeps every junction dihedral in both states, giving
    the single-/dual-anchored redundant-dihedral constraints (cases #2/#3) that
    break factorization.

    Rule (p. 4409): keep the dihedrals reaching into the dummy group from ONE
    physical anchor atom; delete (zero, in the dummy state) every other junction
    dihedral into the group.  The kept anchor dihedrals are left *unchanged* --
    their real periodic terms already anchor the group's rotation with a
    torsional profile identical to the physical atom's, so the group's pose is
    well defined (no flapping) and the physical subsystem is untouched.

    Boresch (Table 3 fn b) additionally *replaces* the kept dihedral with a stiff
    single-minimum term to guarantee no flapping.  We deliberately do NOT do that
    here: the gentop stage has no coordinates, so a hard anchor can only be placed
    at a fixed phase, and empirically (anchor force-constant sweep in scratchpad/)
    a phase that does not match the built geometry couples the physical bridge to
    the dummy and displaces the physical heavy atoms by tens of pm -- a larger
    error than the redundancy it removes.  Keeping the real dihedral avoids that
    while still deleting the factorization-breaking redundant terms.

    Only strictly terminal, single-branch junctions are handled here; dual/triple/
    ring junctions (more than one physical bridge, a group atom bonded to two
    physical atoms, or a bridge carrying a second dummy group) are left as-is
    (they need the fuller Boresch section-2.3 treatment and are over-anchored but
    flapping-free under the naive scheme).
    """
    # bond adjacency over atom ids
    adj = {}
    for b in topol.bonds:
        adj.setdefault(b[0].id, set()).add(b[1].id)
        adj.setdefault(b[1].id, set()).add(b[0].id)
    id2atom = {a.id: a for a in topol.atoms}

    n_del = 0
    n_anchor = 0
    for state, slot in (('A', 5), ('B', 6)):
        dumids = set(a.id for a in topol.atoms if _dummy_in_state(a, state))
        if not dumids:
            continue
        for group in _dc_connected_components(dumids, adj):
            # physical bridge atom(s): physical atoms bonded to the group
            bridges = set()
            for did in group:
                for nb in adj.get(did, ()):
                    if nb not in group and not _dummy_in_state(id2atom[nb], state):
                        bridges.add(nb)
            # terminal junction only: exactly one physical bridge, no ring back
            if len(bridges) != 1:
                continue
            ring = any(
                len([nb for nb in adj.get(did, ())
                     if nb not in group and not _dummy_in_state(id2atom[nb], state)]) > 1
                for did in group)
            if ring:
                continue
            # single dummy branch only: the bridge must not carry a second dummy
            # group (that is a two-branch dual junction -- needs the fuller
            # section-2.3 inter-group treatment, out of scope here).
            X = next(iter(bridges))
            if any(nb not in group for nb in adj.get(X, ())
                   if _dummy_in_state(id2atom[nb], state)):
                continue

            # junction dihedrals reaching into this group in this state
            jdih = []
            for d in topol.dihedrals:
                if len(d) < 7:
                    continue
                params = d[slot]
                if not isinstance(params, list):   # 'NULL' / fake dihedral
                    continue
                ids = [d[0].id, d[1].id, d[2].id, d[3].id]
                in_grp = [i for i in ids if i in group]
                phys = [i for i in ids if not _dummy_in_state(id2atom[i], state)]
                if in_grp and phys:                # spans the junction
                    jdih.append(d)
            if not jdih:
                continue

            # physical "anchor" candidates: the physical END atom (position 0 or
            # 3) of each junction dihedral -- the atom farthest from the dummy.
            def phys_end(d):
                ids = [d[0].id, d[1].id, d[2].id, d[3].id]
                if ids[0] not in group and not _dummy_in_state(id2atom[ids[0]], state):
                    if ids[3] in group or _dummy_in_state(id2atom[ids[3]], state):
                        return ids[0]
                if ids[3] not in group and not _dummy_in_state(id2atom[ids[3]], state):
                    if ids[0] in group or _dummy_in_state(id2atom[ids[0]], state):
                        return ids[3]
                # both ends physical (dummy in the middle) -> use the heavier end
                return None

            # choose the anchor: prefer a heavy (non-H) physical atom that anchors
            # the most junction dihedrals; deterministic tie-break by atom id.
            from collections import Counter
            counts = Counter()
            for d in jdih:
                e = phys_end(d)
                if e is not None:
                    counts[e] += 1
            if not counts:
                continue

            def anchor_key(aid):
                a = id2atom[aid]
                heavy = 0 if (a.atomtype or '').upper().startswith('H') \
                    or (a.name or '').upper().startswith('H') else 1
                return (heavy, counts[aid], -aid)   # heavy > frequent > low-id
            anchor = max(counts, key=anchor_key)

            # Keep the dihedrals reaching into the group from the ONE anchor atom
            # *unchanged* (their real, geometry-consistent periodic terms already
            # anchor the group's rotation without coupling it to a second physical
            # atom); delete every other junction dihedral in the dummy state.
            #
            # We deliberately do NOT replace the kept dihedral with an artificial
            # stiff single-minimum term.  Empirically (anchor force-constant sweep,
            # scratchpad/), a hard anchor at a phase that does not match the built
            # geometry -- which is unavoidable here, since the gentop stage has no
            # coordinates -- couples the physical bridge to the dummy and displaces
            # the physical heavy atoms by tens of pm (worse than the redundancy it
            # removes).  Keeping the real dihedral leaves the group's torsional
            # profile identical to the physical atom's, which is well-defined (no
            # flapping) and does not perturb the physical subsystem.
            anchor_tuples = set()
            for d in jdih:
                func = d[4]
                tup = (d[0].id, d[1].id, d[2].id, d[3].id)
                if phys_end(d) == anchor and func in (1, 4, 9):
                    anchor_tuples.add(tup)                       # keep as-is
                    continue
                d[slot] = _dc_zero_dih(func, d[slot])
                n_del += 1
                if verbose:
                    print('log_> decouple-dummies: deleted junction dihedral '
                          '(state %s) %s-%s-%s-%s' % (state, d[0].name,
                          d[1].name, d[2].name, d[3].name))
            n_anchor += len(anchor_tuples)
    return n_del, n_anchor


def _dc_zero_dih(func, params):
    """Return a zeroed dihedral parameter list for *func* (deletes the term)."""
    if func == 3:                                  # Ryckaert-Bellemans
        return [func, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    if func == 2:                                  # improper (harmonic)
        phi = params[1] if isinstance(params, list) and len(params) > 1 else 0.0
        return [func, phi, 0.0]
    phi = params[1] if isinstance(params, list) and len(params) > 1 else 0.0
    mult = params[-1] if isinstance(params, list) and len(params) > 3 else 1
    return [func, phi, 0.0, mult]                  # proper: zero the force constant


def extract_atoms(residue):
    """Return a name→atom dict for a dual-branch Cys residue (C2CD/D2DC).

    Expects split-atom naming: CB1/CB2 for the two CB branches and SG1/SG2
    for the two SG branches, as defined in the CHARMM36m C2CD/D2DC RTP/MTP.
    CB1/SG1 carry the A-state (free-thiol CYS) parameters;
    CB2/SG2 carry the B-state (disulfide CYS) parameters.
    """
    keys = ('N', 'HN', 'CA', 'HA',
            'CB1', 'CB2', 'HB11', 'HB12', 'HB21', 'HB22',
            'SG1', 'SG2', 'HG1', 'C', 'O')
    atoms = {k: '' for k in keys}
    for a in residue.atoms:
        if a.name in atoms:
            atoms[a.name] = a
    missing = [k for k in ('CA', 'CB1', 'CB2', 'SG1', 'SG2') if atoms[k] == '']
    if missing:
        raise ValueError(
            f"Residue {residue.resname}:{residue.id} missing atoms: {missing}")
    return atoms


def extract_nos_atoms(cys_res, lys_res):
    """Return (cys, lys) name->atom dicts for a NOS hybrid pair.

    Cys side (C2CO, split-atom): CB2/SG2 = B-state branch, OD2 = bridging O.
    Lys side (K2KO, split-atom): CE2/NZ2 = B-state (NOS) branch, HZ21 its proton;
                                 CD (common) for the CD-CE2-NZ2-OD dihedral.
    """
    cwant = ('CA', 'CB2', 'SG2', 'OD2')
    lwant = ('CD', 'CE2', 'NZ2', 'HZ21', 'HE21', 'HE22')
    cys = {k: None for k in cwant}
    lys = {k: None for k in lwant}
    for a in cys_res.atoms:
        if a.name in cys:
            cys[a.name] = a
    for a in lys_res.atoms:
        if a.name in lys:
            lys[a.name] = a
    cmiss = [k for k in ('CB2', 'SG2', 'OD2') if cys[k] is None]
    lmiss = [k for k in ('CE2', 'NZ2') if lys[k] is None]
    if cmiss:
        raise ValueError(f"O2CN residue {cys_res.resname}:{cys_res.id} "
                         f"missing atoms: {cmiss}")
    if lmiss:
        raise ValueError(f"O2KN residue {lys_res.resname}:{lys_res.id} "
                         f"missing atoms: {lmiss}")
    return cys, lys


def _zero_B_if_has_dummy_angle(entry):
    """Zero B-state force constants when any angle atom carries a DUM type."""
    if any('DUM' in a.typeB for a in entry[:3]):
        B = entry[5]
        entry[5] = [B[0], B[1]] + [0.0] * (len(B) - 2)


def _zero_B_if_has_dummy_dihedral(entry):
    """Zero B-state dihedral amplitude when any torsion atom carries a DUM type."""
    if any('DUM' in a.typeB for a in entry[:4]):
        B = entry[6]
        entry[6] = [B[0], B[1], 0.0] + B[3:]


# The bespoke disulfide/NOS bonded-term detectors (_add_ss_bonded_terms /
# _add_nos_terms and their _add_*_bonded_terms drivers) were retired: all
# bond-forming cross-links now flow through _apply_inter_bonds + the graph rule
# induced_bond_morph, driven by the [ inter_bonds ] MTP section or the
# _LEGACY_INTER_BONDS recipe table.  The declarative equivalents disulfide_morph
# / nos_morph live in pmx.covalent.edit.



def _terminal_deletion_decouplings(topol, rlist, verbose=False):
    """Zero B-state force constants for dihedrals spanning the deletion boundary.

    In state B the deletion target is fully dummy. Dihedral terms that span
    the neighbor-deletion boundary are zeroed in state B to remove phantom
    torsional bias on the real atoms.
    """
    for r in rlist:
        if not (r.resname.endswith('deC') or r.resname.endswith('deN')):
            continue

        chain_res = r.chain.residues
        idx = chain_res.index(r)

        if r.resname.endswith('deC') and idx + 1 < len(chain_res):
            del_res = chain_res[idx + 1]
        elif r.resname.endswith('deN') and idx - 1 >= 0:
            del_res = chain_res[idx - 1]
        else:
            continue

        from .covalent.edit import (apply_bond_morph, InterResidueBondMorph,
                                     Interaction)
        neighbor_ids = set(a.id for a in r.atoms)
        del_ids = set(a.id for a in del_res.atoms)

        # collect the dihedrals spanning the deletion boundary, then hand them
        # to the unified applier in 'decouple' mode (zeros their B-state) -- the
        # same zeroing logic as the disulfide/NOS path, one implementation.
        interactions = [
            Interaction('dihedral', [dih[0], dih[1], dih[2], dih[3]])
            for dih in topol.dihedrals
            if ({dih[0].id, dih[1].id, dih[2].id, dih[3].id} & neighbor_ids
                and {dih[0].id, dih[1].id, dih[2].id, dih[3].id} & del_ids)
        ]
        if verbose:
            print('log_> Decoupling cross-residue dihedrals: '
                  '%s (res %d) -- %s (res %d): %d dihedrals'
                  % (r.resname, r.id, del_res.resname, del_res.id,
                     len(interactions)))
        apply_bond_morph(topol, InterResidueBondMorph('decouple', interactions))


def _find_topol_atom(rlist, resid, chain, atom_name):
    """Locate an atom by (resid, chain, name) in the parsed topology residues.

    Prefers the bound-branch split name (``NZ`` -> ``NZ2``) so a user can name the
    plain reacting atom in ``gentop --bond`` without knowing the split naming.
    """
    for r in rlist:
        if r.id != resid:
            continue
        if chain and getattr(r, 'chain_id', '') not in ('', None, chain):
            continue
        by = {a.name: a for a in r.atoms}
        # Try the exact name first (a ligand/partner atom named e.g. "C1"), then
        # the bound-branch split name ("SG" -> "SG2").  Split hybrid residues have
        # no plain reacting-atom name (only X1/X2), so a disulfide/NOS "SG"/"NZ"
        # still falls through to SG2/NZ2; but a covalent-inhibitor ligand atom
        # "C1" resolves to itself instead of grabbing an aromatic "C12".
        for cand in (atom_name, atom_name + '2'):
            if cand in by:
                return by[cand]
    return None


def _apply_inter_bonds(topol, rlist, explicit_bonds, verbose=False):
    """Form the user-declared inter-residue cross-link bonds (``gentop --bond``).

    Each bond is a pair of atoms named by ``(resid, chain, atom)``.  The reacting
    residues must already be clawed back into their dual-topology split (via
    ``pmx mutate --react`` or a shipped split hybrid) so the bond lands on a
    state-A dummy branch and its 1-4 exclusions are harmless.  The induced angles
    and dihedrals follow from the union bond graph -- the same rule GROMACS uses
    to generate bonded terms from bonds -- emitted in dual (A/B) state form
    (:func:`pmx.covalent.edit.induced_bond_morph`).

    This is the single, chemistry-agnostic way an inter-residue bond is added:
    the user names the atom pair, nothing is auto-detected by residue name.
    (A resname->default-bond recognition layer could sit on top of this later.)
    """
    if not explicit_bonds:
        return
    from .covalent.edit import induced_bond_morph, apply_bond_morph

    # union bond graph over the parsed topology (Atom-object adjacency)
    adjacency = {}
    for b in topol.bonds:
        a1, a2 = b[0], b[1]
        adjacency.setdefault(a1, set()).add(a2)
        adjacency.setdefault(a2, set()).add(a1)

    # resolve each --bond pair to its two topology atoms (plain name -> *2 branch).
    # Search ALL residues of this molecule, not just the clawed-back hybrids: a
    # targeted-covalent-inhibitor bond joins a hybrid Cys to a NON-hybrid ligand
    # residue, which is absent from the hybrid-only rlist.  A bond belonging to
    # another chain still has neither endpoint here and is skipped below.
    all_res = getattr(topol, 'residues', None) or rlist
    applied, forming = set(), []
    for eb in explicit_bonds:
        a = _find_topol_atom(all_res, eb['resid_a'], eb.get('chain_a'), eb['atom_a'])
        b = _find_topol_atom(all_res, eb['resid_b'], eb.get('chain_b'), eb['atom_b'])
        if a is None and b is None:
            # Neither endpoint is in this molecule -- the bond belongs to another
            # chain (e.g. the protein in a doublebox; this rlist is the reference
            # peptide).  process_topol runs per-molecule, so skip here; the bond
            # is applied when the molecule that owns both atoms is processed.  A
            # bond never resolved by ANY molecule is caught by gen_hybrid_top's
            # post-processing unresolved-bond check (eb['_applied']).
            continue
        if a is None or b is None:
            raise RuntimeError(
                "gentop --bond: could not resolve %s:%s or %s:%s in the topology "
                "(one endpoint found, the other not -- an inter-residue bond must "
                "lie within a single molecule/chain)"
                % (eb['resid_a'], eb['atom_a'], eb['resid_b'], eb['atom_b']))
        eb['_applied'] = True
        key = frozenset((a.id, b.id))
        if key in applied:
            continue
        applied.add(key)
        forming.append((a, b, eb.get('mode', 'form')))

    # add all forming edges to the graph *before* enumerating, so terms shared
    # between bonds (a SONOS ODa-NZ-ODb cross-angle, both bonds meeting at NZ)
    # are generated.
    for a, b, _ in forming:
        adjacency.setdefault(a, set()).add(b)
        adjacency.setdefault(b, set()).add(a)

    by_mode = {}
    for a, b, mode in forming:
        by_mode.setdefault(mode, []).append((a, b))
    for mode, bonds in by_mode.items():
        morph = induced_bond_morph(bonds, adjacency, mode=mode)
        if verbose:
            print('log_> gentop --bond: %d bond(s) mode=%s -> %d terms'
                  % (len(bonds), mode, len(morph.interactions)))
        apply_bond_morph(topol, morph)


def _find_bonded_entries(topol, verbose=False):
    count = 0
    for b in topol.bonds:
        try:
            a1, a2, func = b
        except: # for the cases where params explicitly listed
            a1, a2, func, otherParam = b
        A, B = _check_case([a1, a2])
        if a1.atomtypeB is not None or a2.atomtypeB is not None:
            count += 1
            error_occured = False
            astate = None
            bstate = None
            if A == 'AA' and B == 'AA':  # we need A and B state
                astate = topol.BondedParams.get_bond_param(a1.type, a2.type)
                bstate = topol.BondedParams.get_bond_param(a1.typeB, a2.typeB)
                b.extend([astate, bstate])
            elif 'D' in A and B == 'AA':
                bstate = topol.BondedParams.get_bond_param(a1.typeB, a2.typeB)
                astate = bstate
                b.extend([astate, bstate])
            elif 'D' in B and A == 'AA':
                astate = topol.BondedParams.get_bond_param(a1.type, a2.type)
                bstate = astate
                b.extend([astate, bstate])
            # catch errors
            elif 'D' in B and 'D' in A:
                print('Error: fake bond '
                      '%s-%s (%s-%s -> %s-%s)' % (a1.name, a2.name,
                                                  a1.atomtype, a2.atomtype,
                                                  a1.atomtypeB, a2.atomtypeB))
                sys.exit(1)
            if astate is None:
                print('Error A: No bond entry found for astate '
                      '%s-%s (%s-%s -> %s-%s)' % (a1.name, a2.name,
                                                  a1.atomtype, a2.atomtype,
                                                  a1.type, a2.type))
                print('Resi = %s' % a1.resname)
                error_occured = True

            if bstate is None:
                print('Error B: No bond entry found for bstate '
                      '%s-%s (%s-%s -> %s-%s)' % (a1.name, a2.name,
                                                  a1.atomtypeB, a2.atomtypeB,
                                                  a1.typeB, a2.typeB))
                error_occured = True
            if error_occured:
                sys.exit(1)
    if verbose is True:
        print('log_> Making bonds for state B -> %d bonds with perturbed atoms'
              % count)


def _find_angle_entries(topol, verbose=False):
    count = 0
    for a in topol.angles:
        try:
            a1, a2, a3, func = a
        except: # for the cases where params explicitly listed
            a1, a2, a3, func, otherParam = a
        astate = None
        bstate = None
        if a1.atomtypeB is not None or \
           a2.atomtypeB is not None or \
           a3.atomtypeB is not None:
            A, B = _check_case([a1, a2, a3])
            count += 1
            if 'D' in A and 'D' in B:  # fake angle
                if func == 5:
                    astate = [1, 0, 0, 0, 0]
                    bstate = [1, 0, 0, 0, 0]
                else:
                    astate = [1, 0, 0]
                    bstate = [1, 0, 0]
                a.extend([astate, bstate])
            elif A == 'AAA' and 'D' in B:   # take astate
                astate = topol.BondedParams.get_angle_param(a1.type,
                                                            a2.type,
                                                            a3.type)
                bstate = astate
                a.extend([astate, bstate])
            elif 'D' in A and B == 'AAA':   # take bstate
                bstate = topol.BondedParams.get_angle_param(a1.typeB,
                                                            a2.typeB,
                                                            a3.typeB)
                astate = bstate
                a.extend([astate, bstate])
            elif A == 'AAA' and B == 'AAA':
                if a1.atomtypeB != a1.atomtype or \
                   a2.atomtypeB != a2.atomtype or \
                   a3.atomtypeB != a3.atomtype:
                    astate = topol.BondedParams.get_angle_param(a1.type,
                                                                a2.type,
                                                                a3.type)
                    bstate = topol.BondedParams.get_angle_param(a1.typeB,
                                                                a2.typeB,
                                                                a3.typeB)
                    a.extend([astate, bstate])
                else:
                    astate = topol.BondedParams.get_angle_param(a1.type,
                                                                a2.type,
                                                                a3.type)
                    bstate = astate
                    a.extend([astate, bstate])
            if astate is None:
                raise MissingTopolParamError("No angle entry (state A)",
                                             [a1, a2, a3])
            if bstate is None:
                raise MissingTopolParamError("No angle entry (state B)",
                                             [a1, a2, a3])

    if verbose is True:
        print('log_> Making angles for state B -> %d angles with perturbed atoms'
              % count)


def _find_dihedral_entries(topol, rlist, rdic, dih_predef_default,
                           verbose=False, scaleDih=1.0):
    count = 0
    nfake = 0
    # here I will accumulate multiple entries of type 9 dihedrals
    dih9 = []
    # need to store already visited angles to avoid multiple re-definitions
    visited_dih = []
    # scale dihedrals with dummies
    #    scaleDih = 0.0

    for d in topol.dihedrals:
        if len(d) >= 6:
            # only consider the dihedral, if it has not been encountered so far
            undef = 0
            encountered = 0

            undef, encountered = _is_dih_undef(dih_predef_default, d)
            encountered = _is_dih_encountered_strict(visited_dih, d,
                                                     encountered)
            if(encountered == 1):
                continue

            visited_dih.append(d)

            if undef == 3 or undef == 4:
                a1 = d[0]
                a2 = d[1]
                a3 = d[2]
                a4 = d[3]
                val = ''
                func = 9
            elif(len(d) == 6):
                a1, a2, a3, a4, func, val = d
            else:
                a1, a2, a3, a4, func, val, rest = d

            backup_d = d[:6]

            if _atoms_morphe([a1, a2, a3, a4]):
                A, B = _check_case(d[:4])
                if A != 'AAAA' and B != 'AAAA':
                    nfake += 1
                    # these are fake dihedrals
                    d[5] = 'NULL'
                    d.append('NULL')
                else:
                    count += 1
                    astate = []
                    bstate = []
                    if A == 'AAAA' and B != 'AAAA':
                        foo = topol.BondedParams.get_dihedral_param(a1.type,
                                                                    a2.type,
                                                                    a3.type,
                                                                    a4.type,
                                                                    func)
                        # need to check if the dihedral has torsion pre-defined
                        counter = _check_dih_ILDN_OPLS(topol, rlist, rdic,
                                                       a1, a2, a3, a4)
                        if counter == 42:
                            continue

                        for ast in foo:
                            if counter == 0:
                                d[5] = ast
                                bst = _scale_dih( ast, scaleDih )
                                d.append(bst)
                            else:
                                alus = backup_d[:]
                                alus[5] = ast
                                bst = _scale_dih( ast, scaleDih )
                                alus.append(bst)
                                dih9.append(alus)
                            counter = 1

                    elif B == 'AAAA' and A != 'AAAA':
                        foo = topol.BondedParams.get_dihedral_param(a1.typeB,
                                                                    a2.typeB,
                                                                    a3.typeB,
                                                                    a4.typeB,
                                                                    func)
                        # need to check if the dihedral has torsion pre-defined
                        counter = _check_dih_ILDN_OPLS(topol, rlist, rdic, a1, a2, a3, a4)

                        if counter == 42:
                            continue

                        for bst in foo:
                            if counter == 0:
                                ast = _scale_dih( bst, scaleDih )
                                d[5] = ast
                                d.append(bst)
                            else:
                                alus = backup_d[:]
                                ast = _scale_dih( bst, scaleDih )
                                alus[5] = ast
                                alus.append(bst)
                                dih9.append(alus)
                            counter = 1

                    elif A == 'AAAA' and B == 'AAAA':
                        # disappear/appear dihedrals, do not morphe
                        if val == '':
                            if(undef != 4):
                                astate = topol.BondedParams.get_dihedral_param(a1.type,
                                                                               a2.type,
                                                                               a3.type,
                                                                               a4.type,
                                                                               func)
                        else:
                            if(undef != 4):
                                astate = topol.BondedParams.get_dihedral_param(a1.type,
                                                                               a2.type,
                                                                               a3.type,
                                                                               a4.type,
                                                                               func)
                        # if types_morphe([a1,a2,a3,a4]):
                        if (undef != 3):
                                bstate = topol.BondedParams.get_dihedral_param(a1.typeB,
                                                                               a2.typeB,
                                                                               a3.typeB,
                                                                               a4.typeB,
                                                                               func)

                        if undef == 1 and astate == []:
                            continue
                        elif undef == 1 and (astate[0][0] == 4 or astate[0][0] == 2):
                            continue
                        elif undef == 2 and bstate == []:
                            continue
                        elif undef == 2 and (bstate[0][0] == 4 or bstate[0][0] == 2):
                            continue

                        # need to check if the dihedral has torsion pre-defined
                        counter = _check_dih_ILDN_OPLS(topol, rlist, rdic,
                                                       a1, a2, a3, a4)

                        # torsion for both states defined, change nothing
                        if counter == 42:
                            continue

                        # A state disappears (when going to B)
                        for ast in astate:
                            if counter == 0:
                                d[5] = ast
                                if ast[0] == 3:
                                    bst = [ast[0], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                                elif ast[0] == 2:
                                    bst = [ast[0], ast[1], 0.0]
                                else:
                                    bst = [ast[0], ast[1], 0.0, ast[-1]]
                                if len(d) == 6:
                                    d.append(bst)
                                else:
                                    d[6] = bst
                                counter = 1
                            elif counter == 1 or counter == 3:
                                alus = backup_d[:]
                                alus[5] = ast
                                if ast[0] == 3:
                                    bst = [ast[0], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                                elif ast[0] == 2:
                                    bst = [ast[0], ast[1], 0.0]
                                else:
                                    bst = [ast[0], ast[1], 0.0, ast[-1]]

                                alus.append(bst)
                                dih9.append(alus)

                        # B state disappears (when going to A)
                        # all must go to the new (appendable) dihedrals
                        for bst in bstate:
                            if counter == 0:
                                if bst[0] == 3:
                                    ast = [bst[0], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                                elif bst[0] == 2:
                                    ast = [bst[0], bst[1], 0.0]
                                else:
                                    ast = [bst[0], bst[1], 0.0, bst[-1]]
                                d[5] = ast
                                if len(d) == 6:
                                    d.append(bst)
                                else:
                                    d[6] = bst
                                counter = 1
                            elif counter == 1 or counter == 2:
                                alus = backup_d[:]
                                if ast[0] == 3:
                                    ast = [bst[0], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                                elif ast[0] == 2:
                                    ast = [bst[0], bst[1], 0.0]
                                else:
                                    ast = [bst[0], bst[1], 0.0, bst[-1]]
                                alus[5] = ast
                                alus.append(bst)
                                dih9.append(alus)


                    if astate is None:
                        print('Error: No dihedral angle found (state A: predefined state B) for:')
                        print('{0} {1} {2} {3}'.format(a1.resname, a2.resname, a3.resname, a4.resname))
                        print('{0} {1} {2} {3}'.format(a1.name, a2.name, a3.name, a4.name, func))
                        print('{0} {1} {2} {3}'.format(a1.atomtype, a2.atomtype, a3.atomtype, a4.atomtype))
                        print('{0} {1} {2} {3}'.format(a1.type, a2.type, a3.type, a4.type))
                        print('{0} {1} {2} {3}'.format(a1.atomtypeB, a2.atomtypeB, a3.atomtypeB, a4.atomtypeB))
                        print('{0} {1} {2} {3}'.format(a1.typeB, a2.typeB, a3.typeB, a4.typeB))
                        print('{0}'.format(astate))
                        print('{0}'.format(d))
                        sys.exit(1)

                    if bstate is None:
                        print('Error: No dihedral angle found (state B: predefined state A) for:')
                        print('{0} {1} {2} {3}'.format(a1.resname, a2.resname, a3.resname, a4.resname))
                        print('{0} {1} {2} {3}'.format(a1.name, a2.name, a3.name, a4.name, func))
                        print('{0} {1} {2} {3}'.format(a1.atomtype, a2.atomtype, a3.atomtype, a4.atomtype))
                        print('{0} {1} {2} {3}'.format(a1.type, a2.type, a3.type, a4.type))
                        print('{0} {1} {2} {3}'.format(a1.atomtypeB, a2.atomtypeB, a3.atomtypeB, a4.atomtypeB))
                        print('{0} {1} {2} {3}'.format(a1.typeB, a2.typeB, a3.typeB, a4.typeB))
                        print('{0}'.format(d))
                        sys.exit(1)

                    # the previous two if sentences will kill the execution
                    # if anything goes wrong. Therefore, I am not afraid to do
                    # d.append() already in the previous steps


    topol.dihedrals.extend(dih9)
    if verbose is True:
        print('log_> Making dihedrals for state B -> %d dihedrals with perturbed atoms' % count)
        print('log_> Removed %d fake dihedrals' % nfake)


def _check_dih_ILDN_OPLS(topol, rlist, rdic, a1, a2, a3, a4):
    counter = 0
    for r in rlist:
        if r.id == a2.resnr:
            dih = rdic[r.resname][3]
            for d in dih:
                if counter == 1:
                    break
                al = []
                for name in d[:4]:
                    atom = r.fetch(name)[0]
                    al.append(atom)

                if (a1.id == al[0].id and a2.id == al[1].id and
                   a3.id == al[2].id and a4.id == al[3].id) or \
                   (a1.id == al[3].id and a2.id == al[2].id and
                   a3.id == al[1].id and a4.id == al[0].id):
                    # ---------------
                    # checks for OPLS
                    # ---------------
                    if d[4].startswith('dih_') and d[5].startswith('dih_'):
                        counter = 42
                        break
                    elif d[4].startswith('dih_') and 'un' in d[5]:
                        counter = 2
                        break
                    elif 'un' in d[4] and d[5].startswith('dih_'):
                        counter = 3
                        break
                    # ---------------
                    # checks for ILDN
                    # ---------------
                    if d[4].startswith('torsion') and d[5].startswith('torsion'):
                        counter = 42
                        break
                    elif d[4].startswith('torsion') and 'un' in d[5]:
                        counter = 2
                        break
                    elif 'un' in d[4] and d[5].startswith('torsion'):
                        counter = 3
                        break
                    elif 'un' in d[4]:
                        counter = 1
                        break
    return counter


def _is_dih_undef(visited_dih, d):
    encountered = 0
    undef = 0
    for dih in visited_dih:
        if (dih[0].id == d[0].id and dih[1].id == d[1].id
           and dih[2].id == d[2].id and dih[3].id == d[3].id):
            if (dih[-1] == 'undefA'):
                undef = 1
                break
            elif (dih[-1] == 'undefB'):
                undef = 2
                break
            elif (dih[-1] == 'undefA_ildn'):
                undef = 3
                break
            elif (dih[-1] == 'undefB_ildn'):
                undef = 4
                break
            else:
                encountered = 1
                break

    return (undef, encountered)


def _is_dih_encountered_strict(visited_dih, d, encountered):
    for dih in visited_dih:
        if (dih[0].id == d[0].id and dih[1].id == d[1].id
           and dih[2].id == d[2].id and dih[3].id == d[3].id):
            encountered = 1
            break
    return encountered


def _get_torsion_multiplicity(name):
    foo = list(name)
    mult = int(foo[-1])
    return mult


def _explicit_defined_dihedrals(filename, ff):
    """ get the #define dihedral entries explicitly for ILDN
    """
    lines = open(filename).readlines()
    output = {}
    ffnamelower = ff.lower()
    for line in lines:
        if line.startswith('#define'):
            entr = line.split()
            name = entr[1]
            if "amber" in ffnamelower: #ffnamelower.startswith('amber'):
                params = [9, float(entr[2]), float(entr[3]), int(entr[4])]
            elif ffnamelower.startswith('opls'):
                if len(entr) == 8:  # dihedral
                    params = [3, float(entr[2]), float(entr[3]), float(entr[4]),
                              float(entr[5]), float(entr[6]), float(entr[7])]
                elif len(entr) == 5:
                    params = [1, float(entr[2]), float(entr[3]), float(entr[4])]
            output[name] = params
    return output


def _is_ildn_dih_encountered(ildn_used, d, encountered):
    for dih in ildn_used:
        if (dih[0] == d[0] and dih[1] == d[1] and dih[2] == d[2]
           and dih[3] == d[3] and dih[4] == d[4]):
            encountered = 1
    return encountered


def _find_predefined_dihedrals(topol, rlist, rdic, ffbonded,
                               dih_predef_default, ff, scaleDih=1.0):

    dih9 = []  # here I will accumulate multiple entries of type 9 dihedrals
    explicit_def = _explicit_defined_dihedrals(ffbonded, ff)
    ildn_used = []  # ildn dihedrals that already were encountered
    opls_used = []  # opls dihedrals that already were encountered
#    scaleDih = 0.0  # scale dihedrals with dummies

    for r in rlist:
        dih = rdic[r.resname][3]
        imp = rdic[r.resname][2]
        for d in imp+dih:
            al = []
            for name in d[:4]:
                if name.startswith('+'):
                    next = r.chain.fetch_residue(idx=r.id+1)
#                    next = r.chain.residues[idx+1]
                    atom = next.fetch(name[1:])[0]
                    al.append(atom)
                elif name.startswith('-'):
                    prev = r.chain.fetch_residue(idx=r.id-1)
#                    prev = r.chain.residues[idx-1]
                    atom = prev.fetch(name[1:])[0]
                    al.append(atom)
                else:
                    atom = r.fetch(name)[0]
                    al.append(atom)

            for dx in topol.dihedrals:
                func = dx[4]
                if (dx[0].id == al[0].id and
                   dx[1].id == al[1].id and
                   dx[2].id == al[2].id and
                   dx[3].id == al[3].id) or \
                   (dx[0].id == al[3].id and
                   dx[1].id == al[2].id and
                   dx[2].id == al[1].id and
                   dx[3].id == al[0].id):
                    # the following checks are needed for amber99sb*-ildn
                    # do not overwrite proper (type9) with improper (type4)
                    if 'default-star' in d[4] and dx[4] == 9:
                        continue
                    # is the dihedral already found for ILDN
                    encountered = 0
                    foobar = [dx[0].id, dx[1].id, dx[2].id, dx[3].id, d[4]]
                    encountered = _is_ildn_dih_encountered(ildn_used, foobar, encountered)
                    foobar = [dx[0].id, dx[1].id, dx[2].id, dx[3].id, d[5]]
                    encountered = _is_ildn_dih_encountered(ildn_used, foobar, encountered)
                    if encountered == 1:
                        continue
                    if 'tors' in d[4]:
                        if 'tors' not in dx[5]:
                            continue

                    # check for opls
                    encountered = 0
                    foobar = [dx[0].id, dx[1].id, dx[2].id, dx[3].id, d[4]]
                    encountered = _is_ildn_dih_encountered(opls_used, foobar, encountered)
                    foobar = [dx[0].id, dx[1].id, dx[2].id, dx[3].id, d[5]]
                    encountered = _is_ildn_dih_encountered(opls_used, foobar, encountered)
                    if encountered == 1:
                        continue
                    if 'dih_' in d[4]:
                        if 'dih_' not in dx[5]:
                            continue

                    A, B = _check_case(al[:4])
                    paramA = topol.BondedParams.get_dihedral_param(al[0].type,
                                                                   al[1].type,
                                                                   al[2].type,
                                                                   al[3].type,
                                                                   func)
                    paramB = topol.BondedParams.get_dihedral_param(al[0].typeB,
                                                                   al[1].typeB,
                                                                   al[2].typeB,
                                                                   al[3].typeB,
                                                                   func)

                    astate = []
                    bstate = []
                    backup_dx = dx[:]
                    backup_dx2 = dx[:]

                    if d[4] == 'default-A':  # amber99sb
                        if 'un' in d[5]:
                            backup_dx2.append('undefB')
                        dih_predef_default.append(backup_dx2)
                        astate = paramA
                    elif d[4] == 'default-B':  # amber99sb
                        if 'un' in d[5]:
                            backup_dx2.append('undefB')
                        dih_predef_default.append(backup_dx2)
                        astate = paramB
                    elif d[4] == 'default-star':  # amber99sb*
                        foo = [4, 105.4, 0.75, 1]
                        astate.append(foo)
                    elif d[4].startswith('torsion_'):  # amber99sb-ildn
                        if 'un' in d[5]:
                            backup_dx2.append('undefB_ildn')
                        dih_predef_default.append(backup_dx2)
                        foo = explicit_def[d[4]]
                        astate.append(foo)
                        foobar = [dx[0].id, dx[1].id, dx[2].id, dx[3].id, d[4]]
                        ildn_used.append(foobar)
                        func = 9
                    elif d[4].startswith('dih_'):  # opls proper
                        if 'un' in d[5]:
                            backup_dx2.append('undefB')
                        dih_predef_default.append(backup_dx2)
                        foo = explicit_def[d[4]]
                        astate.append(foo)
                        foobar = [dx[0].id, dx[1].id, dx[2].id, dx[3].id, d[4]]
                        opls_used.append(foobar)
                        func = 3
                    elif d[4].startswith('improper_'):  # opls improper
                        if 'un' in d[5]:
                            backup_dx2.append('undefB')
                        dih_predef_default.append(backup_dx2)
                        foo = explicit_def[d[4]]
                        astate.append(foo)
                        foobar = [dx[0].id, dx[1].id, dx[2].id, dx[3].id, d[4]]
                        opls_used.append(foobar)
                        func = 1
                    elif 'un' in d[4]:  # amber99sb-ildn and opls and others
                        if d[5].startswith('torsion_'):
                            backup_dx2.append('undefA_ildn')
                        else:
                            backup_dx2.append('undefA')
                        dih_predef_default.append(backup_dx2)
                        astate = ''
                    else:
                        astate = d[4]

                    if d[5] == 'default-A':  # amber99sb
                        bstate = paramA
                    elif d[5] == 'default-B':  # amber99sb
                        bstate = paramB
                    elif d[5] == 'default-star':  # amber99sb*
                        foo = [4, 105.4, 0.75, 1]
                        bstate.append(foo)
                    elif d[5].startswith('torsion_'):  # amber99sb-ildn
                        foo = explicit_def[d[5]]
                        bstate.append(foo)
                        foobar = [dx[0].id, dx[1].id, dx[2].id, dx[3].id, d[5]]
                        ildn_used.append(foobar)
                        func = 9
                    elif d[5].startswith('dih_'):  # opls proper
                        foo = explicit_def[d[5]]
                        bstate.append(foo)
                        foobar = [dx[0].id, dx[1].id, dx[2].id, dx[3].id, d[5]]
                        opls_used.append(foobar)
                        func = 3
                    elif d[5].startswith('improper_'):  # opls improper
                        foo = explicit_def[d[5]]
                        bstate.append(foo)
                        foobar = [dx[0].id, dx[1].id, dx[2].id, dx[3].id, d[5]]
                        ildn_used.append(foobar)
                        func = 1
                    elif 'un' in d[5]:  # amber99sb-ildn and opls and others
                        bstate = ''
                    else:
                        bstate = d[5]

                    # A state
                    counter = 0
                    for foo in astate:
                        if counter == 0:
                            dx[4] = func
                            if 'D' in A:
                               foo =  _scale_dih( foo, scaleDih )
                            dx[5] = foo
                            if(func == 3):
                                bar = [foo[0], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                            elif(func == 2):
                                bar = [foo[0], foo[1], 0.0]
                            else:
                                bar = [foo[0], foo[1], 0.0, foo[-1]]
                            dx.append(bar)
                        else:
                            alus = backup_dx[:]
                            if 'D' in A:
                               foo =  _scale_dih( foo, scaleDih )
                            alus[5] = foo
                            if func == 3:
                                bar = [foo[0], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                            elif(func == 2):
                                bar = [foo[0], foo[1], 0.0]
                            else:
                                bar = [foo[0], foo[1], 0.0, foo[-1]]
                            alus.append(bar)
                            dih9.append(alus)
                        counter = 1

                    # B state
                    for foo in bstate:
                        if counter == 0:
                            dx[4] = func
                            if(func == 3):
                                bar = [foo[0], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                            elif(func == 2):
                                bar = [foo[0], foo[1], 0.0]
                            else:
                                bar = [foo[0], foo[1], 0.0, foo[-1]]
                            dx[5] = bar
                            if 'D' in B:
                               foo =  _scale_dih( foo, scaleDih )
                            dx.append(foo)
                        else:
                            alus = backup_dx[:]
                            if 'D' in B:
                               foo =  _scale_dih( foo, scaleDih )
                            if(func == 3):
                                bar = [foo[0], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                            elif(func == 2):
                                bar = [foo[0], foo[1], 0.0]
                            else:
                                bar = [foo[0], foo[1], 0.0, foo[-1]]
                            alus[5] = bar
                            alus.append(foo)
                            dih9.append(alus)
                        counter = 1

    topol.dihedrals.extend(dih9)


def _get_hybrid_residue(residue_name, mtp_file='ffamber99sb.mtp',
                        version='new', verbose=False):
    if verbose is True:
        print('log_> Scanning database for %s ' % residue_name)
    resi, bonds, imps, diheds, rotdic = read_mtp_entry(residue_name,
                                                       filename=mtp_file,
                                                       version=version)
    if len(resi.atoms) == 0:
        raise mtpError("Hybrid residue %s not found in %s"
                       % (residue_name, mtp_file))
    return resi, bonds, imps, diheds, rotdic


def _inject_supplement_bonded(pmxtop, supplement_bonded_files):
    """Prepend bonded parameters from supplement ITP files into *pmxtop*.

    Supplement ITP files are written by :func:`~pmx.nsaa.prepare_nsaa_ff`
    and contain ``[ bondtypes ]``, ``[ angletypes ]``, and
    ``[ dihedraltypes ]`` sections (plus optionally ``[ atomtypes ]``) for
    AMBER↔GAFF cross-boundary terms and ITP-internal GAFF bonded parameters.

    Entries are prepended to ``pmxtop.BondedParams.bondtypes`` etc. so they
    take priority over the built-in force-field database when
    :func:`_find_bonded_entries`, :func:`_find_angle_entries`, and
    :func:`_find_predefined_dihedrals` perform their lookups.

    Parameters
    ----------
    pmxtop : Topology
    supplement_bonded_files : list of str
    """
    from .nsaa import (_parse_ffbonded_section, _parse_atomtypes_section)

    def _to_typed(entry, n_types):
        """Convert a _parse_ffbonded_section dict to a BondedParser list."""
        types  = entry['types']     # list of str, len == n_types
        funct  = entry['funct']     # int
        params = entry['params']    # list of str
        # BondedParser stores entries as mixed lists: [str*n_types, int, float*n]
        try:
            float_params = [float(p) for p in params]
        except ValueError:
            float_params = params   # leave as-is if non-numeric (shouldn't happen)
        return types + [funct] + float_params

    for supp_file in supplement_bonded_files:
        if not os.path.isfile(supp_file):
            print('log_> WARNING: supplement bonded file not found: %s' % supp_file)
            continue

        # --- [ atomtypes ] ---
        try:
            at_entries = _parse_atomtypes_section(supp_file)
        except Exception:
            at_entries = []
        for e in at_entries:
            name = e[0]
            if name not in pmxtop.NBParams.atomtypes:
                pmxtop.NBParams.atomtypes[name] = {
                    'bond_type': name,
                    'mass':  float(e[2]),
                    'sigma': float(e[5]),
                    'eps':   float(e[6]),
                }

        # --- [ bondtypes ] ---
        bond_entries = _parse_ffbonded_section(supp_file, 'bondtypes', 2)
        for e in reversed(bond_entries):
            pmxtop.BondedParams.bondtypes.insert(0, _to_typed(e, 2))

        # --- [ angletypes ] ---
        angle_entries = _parse_ffbonded_section(supp_file, 'angletypes', 3)
        for e in reversed(angle_entries):
            pmxtop.BondedParams.angletypes.insert(0, _to_typed(e, 3))

        # --- [ dihedraltypes ] ---
        dih_entries = _parse_ffbonded_section(supp_file, 'dihedraltypes', 4)
        for e in reversed(dih_entries):
            pmxtop.BondedParams.dihedraltypes.insert(0, _to_typed(e, 4))

        n_b = len(bond_entries)
        n_a = len(angle_entries)
        n_d = len(dih_entries)
        print('log_> Supplement bonded: %d bonds, %d angles, %d dihedrals '
              'injected from %s' % (n_b, n_a, n_d, supp_file))


def _read_mtp_names(filename):
    """Return all top-level residue names defined in an MTP file.

    Scans for ``[ NAME ]`` section headers and filters out the known
    sub-section keywords (atoms, bonds, morphes, …).

    Parameters
    ----------
    filename : str
        Path to an MTP file.

    Returns
    -------
    names : list of str
        Residue names found in *filename*, in order of appearance.
    """
    _sub_sections = frozenset((
        'morphes', 'atoms', 'bonds', 'impropers',
        'dihedrals', 'rotations', 'coords',
    ))
    names = []
    try:
        with open(filename) as fh:
            for line in fh:
                stripped = line.strip()
                if stripped.startswith('['):
                    # strip inline comments (;) and brackets
                    core = stripped.split(';')[0].strip()
                    if core.startswith('[') and core.endswith(']'):
                        name = core[1:-1].strip()
                        if name and name not in _sub_sections:
                            names.append(name)
    except OSError as exc:
        raise OSError("Cannot read MTP file '%s': %s" % (filename, exc)) from exc
    return names


def _get_hybrid_residues(m, ff, version='new', verbose=False,
                         extra_mtp_files=None):
    """Collect hybrid residues from *m* and load their MTP parameters.

    Parameters
    ----------
    m : Model
    ff : str
        Force-field name.
    version : str
        MTP format version ('new' or 'old').
    verbose : bool
    extra_mtp_files : list of str, optional
        Additional MTP files to search for residue parameters.  These take
        precedence over the force-field's built-in MTP database.  Use this
        to supply NSAA hybrid MTP files generated by :func:`~pmx.nsaa.mutate_nsaa`.

    Returns
    -------
    rlist : list of Residue
    rdic  : dict  {resname: (mol, bonds, imps, diheds, rotdic)}
    """
    # Build resname → mtp_path lookup from extra files (first hit wins)
    extra_name_map = {}
    if extra_mtp_files:
        for mtp_path in extra_mtp_files:
            for name in _read_mtp_names(mtp_path):
                if name not in extra_name_map:
                    extra_name_map[name] = mtp_path

    # Auto-discover the force field's nsaa.mtp: the generative engine
    # (covalent.mutate) accumulates every runtime-built hybrid there, so
    # gen_hybrid_top resolves them with no --extra_mtp.  Explicit extra_mtp_files
    # still win (added first, first-hit-wins).
    from .utils import get_ff_path
    nsaa_mtp = os.path.join(get_ff_path(ff), 'nsaa.mtp')
    if os.path.isfile(nsaa_mtp):
        for name in _read_mtp_names(nsaa_mtp):
            extra_name_map.setdefault(name, nsaa_mtp)

    rdic = {}
    rlist = []
    for res in m.residues:
        if not res.is_hybrid() and res.resname not in extra_name_map:
            continue
        rlist.append(res)
        # Extra MTP files take precedence over the FF database
        if res.resname in extra_name_map:
            mtp_file = extra_name_map[res.resname]
        else:
            mtp_file = get_mtp_file(res, ff)
        mtp = _get_hybrid_residue(residue_name=res.resname,
                                  mtp_file=mtp_file,
                                  version=version,
                                  verbose=verbose)
        rdic[res.resname] = mtp
        hybrid_res = mtp[0]
        # Match each mtp atom to the structure atom BY NAME.  fetchm() silently
        # skips names it cannot find, which used to misalign a positional index
        # (atoms[i]) and raise an opaque IndexError.  Matching by name instead
        # gives a clear, actionable error when a required atom is genuinely absent
        # (e.g. the hybrid was built from a template missing that atom).
        by_name = {a.name: a for a in res.atoms}
        for atom in hybrid_res.atoms:
            target = by_name.get(atom.name)
            if target is None:
                raise mtpError(
                    "hybrid residue %s%s: atom '%s' listed in the mtp is missing "
                    "from the structure topology. The hybrid was likely built from "
                    "a template that lacks this atom (e.g. a crystal template with "
                    "no hydrogens). Rebuild the target residue with a complete-atom "
                    "template." % (res.resname, res.id, atom.name))
            target.atomtypeB = atom.atomtypeB
            target.qB = atom.qB
            target.mB = atom.mB
    return rlist, rdic


def _add_extra_DNA_RNA_impropers(topol, rlist, func_type, stateA, stateB):
    extra_impropers = []
    for r in rlist:
        if r.resname in ['DAT', 'DAC', 'DGC', 'DGT', 'RAU', 'RAC', 'RGC', 'RGU']:
            print('Adding extra improper dihedrals for residue %d-%s'
                  % (r.id, r.resname))
            alist = r.fetchm(['C1\'', 'N9', 'C8', 'DC2'])
            extra_impropers.append(alist)
            alist = r.fetchm(['C1\'', 'N9', 'C4', 'DC6'])
            extra_impropers.append(alist)
        elif r.resname in ['DTA', 'DTG', 'DCG', 'DCA', 'RUA', 'RUG', 'RCG', 'RCA']:
            print('Adding extra improper dihedrals for residue %d-%s'
                  % (r.id, r.resname))
            alist = r.fetchm(['C1\'', 'N1', 'C6', 'DC4'])
            extra_impropers.append(alist)
            alist = r.fetchm(['C1\'', 'N1', 'C2', 'DC8'])
            extra_impropers.append(alist)
    for imp in extra_impropers:
        imp.extend([func_type, [func_type]+stateA, [func_type]+stateB])
    topol.dihedrals += extra_impropers


# =======
# Classes
# =======
class AbsRestraints:
    '''Identifies a set of restraints as defined by Boresch between the
    protein/host and ligand/guest molecule. [5]_ ::

                (Pro)  p3 -- p2          l2   (Lig)
                              \         /  \\
                               p1 --- l1    l3

    The process of choosing atoms involved in the restraints is
    stochastic. Thus, every time you call this function, you might get a
    different set of restraints. If you would like to have a deterministic
    behaviour provide a ``seed``.

    If you have a structure of the protein-ligand complex and want
    to choose restraints based on this, you can do so by providing a
    Model of the ``complex`` as well as the residue name of the ligand
    (``ligname``).

    If you know which atoms to use for the restraints, you can provide their
    atom indices as they are in the input Model objects with the arguments
    ``pro_ids`` and ``lig_ids``. The expected order of atoms is ``l1, l2, l3`` as
    in the scheme above for the ligand, and ``p1, p2, p3`` for the protein.

    Parameters
    ----------
    protein : Model, optional
        model of the protein/host molecule.
    ligand : Model, optional
        model of the ligand/guest molecule.
    complex: Model, optional
        model of the protein-ligand or host-guest complex.
    ligname : str, optional
        name of ligand/host residue. This is required if a complex is provided.
    pro_ids : list, optional
        list with three indeces.
    lig_ids : list, optional
        list with three indeces.
    kbond : float, optional
        force constant to use for the distance restraint.
        Default is 4184 kJ/(mol * nm^2).
    kangle : float, optional
        force constant to use for the angle restraints.
        Default is 41.84 kJ/(mol * rad^2).
    kdihedral : float, optional
        force constant to use for the dihedral restraints.
        Default is 41.84 kJ/(mol * rad^2).
    T : float, optional
        temperature in Kelvin. Default is 298.15 K.
    seed : bool, optional
        random seed.

    Attributes
    ----------
    lig_atoms : list
        list of three ligand/guest Atoms selected.
    pro_atoms : list
        list of three protein/host Atoms selected.
    dist : float
        distance ``l1-p1`` in nm.
    angle1 : float
        angle ``l2-l1-p1`` in degrees.
    angle2 : float
        angle ``l1-p1-p2`` in degrees.
    dihedral1 : float
        dihedral ``l3-l2-l1-p1`` in degrees.
    dihedral2 : float
        dihedral ``l2-l1-p1-p2`` in degrees.
    dihedral3 : float
        dihedral ``l1-p1-p2-p3`` in degrees.
    dg : float
        free energy of restraining the ligand while decoupled. This is calculate
        with equation 32 in Boresch et al. [5]_

    '''
    def __init__(self, protein=None, ligand=None, complex=None, ligname=None,
                 pro_ids=None, lig_ids=None,
                 kbond=4184, kangle=41.84, kdihedral=41.84,
                 T=298.15, seed=None):

        # check which scenario we have: either protein/ligand models are
        # provided separately, or a complex is provided along with the
        # residue name for the ligand.
        self.complex = deepcopy(complex)
        if self.complex is not None:
            if protein is not None or ligand is not None:
                raise ValueError('both arguments "complex" and "protein" or '
                                 '"ligand" have been provided. You should provide '
                                 'either separate protein and ligand Models '
                                 'or a single Model of the complex, but not '
                                 'both')
            if ligname is None:
                raise ValueError('when providing a protein-ligand complex '
                                 'you also need to provide the residue name '
                                 'of the ligand')
            self.ligname = ligname
            self.protein, self.ligand = self._split_complex_model()
        if self.complex is None:
            self.protein = deepcopy(protein)
            self.ligand = deepcopy(ligand)

        # other input parameters
        self.pro_ids = pro_ids
        self.lig_ids = lig_ids
        self.kbond = float(kbond)
        self.kangle = float(kangle)
        self.kdihedral = float(kdihedral)
        self.T = float(T)
        self.seed = seed
        # restraints data
        self.lig_atoms = []
        self.pro_atoms = []
        self.dist = None
        self.angle1 = None
        self.angle2 = None
        self.dihedral1 = None
        self.dihedral2 = None
        self.dihedral3 = None

        # check protein and ligand model are using the same units
        assert self.ligand.unity == self.protein.unity
        if self.ligand.unity == 'A':
            self.f = 10.0
        elif self.ligand.unity == 'nm':
            self.f = 1.0

        # pick atoms and define restraints
        self.select_restraints()

    def select_restraints(self):
        '''Automated restraints selection. For the ligand/guest, three heavy
        atoms close to its center of mass are selected. For the protein,
        the backbone N-CA-C atoms closest to the selected ligand atoms are
        picked. If the protein/host is another small molecule, then three
        atoms close to the ligand atoms are selected.
        '''
        # ===================
        # Choose Ligand Atoms
        # ===================

        # pick 3 atoms near the centre of mass
        if self.lig_ids is None:
            # find centre of mass
            com = Atom()
            com.x = self.ligand.com(vector_only=True)
            # calc distances from com, but exclude H atoms from selection
            distances = [(a, a-com) for a in self.ligand.atoms if a.symbol != 'H']
            # first atom picked is the one closest to centre of mass
            l1 = min(distances, key=lambda x: x[1])[0]
            # make a list of atoms within 3A of atom 1, excluding hydrogens
            atoms = [a for a in self.ligand.atoms if l1.id != a.id and l1-a < 0.3*self.f and a.symbol != 'H']
            # then pick two of them at random
            if self.seed is not None:
                np.random.seed(self.seed)
            l2, l3 = np.random.choice(atoms, size=2)
        else:
            if len(self.lig_ids) != 3:
                raise ValueError()
            l1 = self.ligand.fetch_atoms(key=self.lig_ids[0], how='byid')[0]
            l2 = self.ligand.fetch_atoms(key=self.lig_ids[1], how='byid')[0]
            l3 = self.ligand.fetch_atoms(key=self.lig_ids[2], how='byid')[0]

        self.lig_atoms = [l1, l2, l3]

        # ====================
        # Choose Protein Atoms
        # ====================
        if self.pro_ids is None:
            # check there are atoms close to the ligand (within 1 nm)
            atoms = [a for a in self.protein.atoms if a.symbol != 'H' and l1-a < 1.0*self.f]
            assert len(atoms) > 2
            # if host is a protein, choose backbone atoms
            # -------------------------------------------
            if self.protein.moltype == 'protein':
                N_atoms = [(a, l1-a) for a in atoms if a.name == 'N']
                # get closest N atom
                p1 = min(N_atoms, key=lambda x: x[1])[0]
                # get CA and C of the same residue
                mol = p1.molecule
                p2 = mol.fetch_atoms(key='CA', how='byname')[0]
                p3 = mol.fetch_atoms(key='C', how='byname')[0]

            # if host is not a protein, choose closest atom, plus 2 other random
            # ------------------------------------------------------------------
            else:
                atoms = [(a, l1-a) for a in self.protein.atoms if a.symbol != 'H']
                # first atom is the one closest to l1
                p1 = min(atoms, key=lambda x: x[1])[0]
                # make a list of atoms within 3A of atom 1, excluding hydrogens
                atoms = [a for a in self.protein.atoms if p1.id != a.id and p1-a < 0.3*self.f and a.symbol != 'H']
                # then pick two of them at random
                if self.seed is not None:
                    np.random.seed(self.seed)
                p2, p3 = np.random.choice(atoms, size=2)
        else:
            if len(self.pro_ids) != 3:
                raise ValueError()
            p1 = self.protein.fetch_atoms(key=self.pro_ids[0], how='byid')[0]
            p2 = self.protein.fetch_atoms(key=self.pro_ids[1], how='byid')[0]
            p3 = self.protein.fetch_atoms(key=self.pro_ids[2], how='byid')[0]

        self.pro_atoms = [p1, p2, p3]

        # --------------------------
        # Calculate distances/angles
        # --------------------------
        self.dist = (l1-p1)/self.f  # make sure it is in nm
        self.angle1 = l1.angle(l2, p1, degree=True)
        self.angle2 = p1.angle(l1, p2, degree=True)
        self.dihedral1 = l3.dihedral(l2, l1, p1, degree=True)
        self.dihedral2 = l2.dihedral(l1, p1, p2, degree=True)
        self.dihedral3 = l1.dihedral(p1, p2, p3, degree=True)

        # calc restraints dg
        self.calc_dg()

    def calc_dg(self, release=False):
        '''Calculates the free energy contribution of the
        restraints given the distance, angles, and dihedrals and their
        force constants as defined in the instance.

        Parameters
        ----------
        release : bool, optional
            whether to calculate the free energy of adding (False) or releasing
            (True) the restraints.
        '''
        V0 = 1.66  # standard volume in nm^3
        R = 8.314472*0.001  # Gas constant (kJ/mol/K)
        RT = R * self.T
        prefactor = (8.0*np.power(np.pi, 2.0) * V0 / (np.power(self.dist, 2.0) *
                     np.sin(self.angle1*np.pi/180.0) * np.sin(self.angle2*np.pi/180.0)))
        fconstants = (np.sqrt(self.kbond*self.kangle*self.kangle*self.kdihedral*self.kdihedral*self.kdihedral) /
                      np.power(2.0*np.pi*RT, 3.0))

        if release is True:
            self.dg = -RT * np.log(prefactor * fconstants)
        elif release is False:
            self.dg = RT * np.log(prefactor * fconstants)

    def make_ii(self, switch_on=True, ligand_first=True):
        '''Returns the data needed to add an intermolecular interactions
        section in a topology file to be used as restraints.

        Parameters
        ----------
        switch_on : bool, optional
            whether you want to switch the restraints on (restraints present in
            state B) or off (restraints present in state A). Default is
            True (switching on).
        ligand_first : bool, optional
            whether in the input gro/pdb file of the protein-ligand complex
            you are preparing the ligand atoms will come before (True) or after
            (False) the protein atoms. This is needed in order to adjust the
            atoms indices correctly. If you provided a Model of the complex
            (self.is_complex is True), this argument is disregarded as indices
            as provided will be used.

        Returns
        -------
        ii : dict
            dictionary with "bonds", "angles", and "dihedrals" keys with
            intermolecular_interactions information that can
            be passed to a Topology object.
        '''

        ii = {}
        ii['bonds'] = []
        ii['angles'] = []
        ii['dihedrals'] = []

        # ----------------
        # Sort out indices
        # ----------------
        l1 = self.lig_atoms[0].id
        l2 = self.lig_atoms[1].id
        l3 = self.lig_atoms[2].id
        p1 = self.pro_atoms[0].id
        p2 = self.pro_atoms[1].id
        p3 = self.pro_atoms[2].id
        if self.complex is None:
            if ligand_first is True:
                nshift = len(self.ligand.atoms)
                p1 += nshift
                p2 += nshift
                p3 += nshift
            elif ligand_first is False:
                nshift = len(self.protein.atoms)
                l1 += nshift
                l2 += nshift
                l3 += nshift

        # -----
        # Bonds
        # -----
        if switch_on is True:
            kA = 0.0
            kB = self.kbond
        elif switch_on is False:
            kA = self.kbond
            kB = 0.0

        ii['bonds'].append([l1, p1, 6, [self.dist, kA, self.dist, kB]])

        # ------
        # Angles
        # ------
        if switch_on is True:
            kA = 0.0
            kB = self.kangle
        else:
            kA = self.kangle
            kB = 0.0

        ii['angles'].append([l2, l1, p1, 1, [self.angle1, kA, self.angle1, kB]])
        ii['angles'].append([l1, p1, p2, 1, [self.angle2, kA, self.angle2, kB]])

        # ---------
        # Dihedrals
        # ---------
        if switch_on is True:
            kA = 0.0
            kB = self.kdihedral
        else:
            kA = self.kdihedral
            kB = 0.0

        ii['dihedrals'].append([l3, l2, l1, p1, 2, [self.dihedral1, kA, self.dihedral1, kB]])
        ii['dihedrals'].append([l2, l1, p1, p2, 2, [self.dihedral2, kA, self.dihedral2, kB]])
        ii['dihedrals'].append([l1, p1, p2, p3, 2, [self.dihedral3, kA, self.dihedral3, kB]])

        return ii

    def write_summary(self, fname='restraints.info'):
        '''Write restraints information, including free energy of restraining,
        to file.

        Parameters
        ----------
        fname : str, optional
            name of output file. Default is "restraints.info".
        '''

        summary = """================== Restraints ==================
atoms (lig)   = {0:>6} {1:>6} {2:>6}
atoms (pro)   = {3:>6} {4:>6} {5:>6}
""".format(self.lig_atoms[0].id, self.lig_atoms[1].id, self.lig_atoms[2].id,
           self.pro_atoms[0].id, self.pro_atoms[1].id, self.pro_atoms[2].id)

        summary += """------------------------------------------------
T             = {T:>14.2f} K
kbond         = {kbond:>14.2f} kJ/(mol * nm^2)
kangles       = {kangle:>14.2f} kJ/(mol * rad^2)
kdihedrals    = {kdihedral:>14.2f} kJ/(mol * rad^2)
r0            = {dist:>14.6f} nm
thA           = {angle1:>14.6f} deg
thB           = {angle2:>14.6f} deg
phiA          = {dihedral1:>14.6f} deg
phiB          = {dihedral2:>14.6f} deg
phiC          = {dihedral3:>14.6f} deg
================================================
dG Restraints = {dg:>14.2f} kJ/mol""".format(**self.__dict__)

        summary += "\n              = {0:>14.2f} kcal/mol\n".format(self.dg/4.184)

        with open(fname, 'w') as f:
            f.write(summary)

    def _split_complex_model(self):
        complex = deepcopy(self.complex)
        # find ligand
        ligands = complex.fetch_residues(key=self.ligname)
        # do a couple of checks
        if len(ligands) == 0:
            raise ValueError('no ligand with residue name %s found in Model '
                             'from file %s' % (self.ligname, self.complex.filename))
        elif len(ligands) > 1:
            raise ValueError('multiple ligands with residue name %s found in Model '
                             'from file %s' % (self.ligname, self.complex.filename))

        # make new Model for ligand
        ligand_mol = ligands[0]
        ligand = Model()
        for a in ligand_mol.atoms:
            ligand.atoms.append(a)
        ligand.unity = ligand_mol.atoms[0].unity

        # make protein by removing ligand residue
        complex.remove_residue(ligands[0], renumber_atoms=False, renumber_residues=False)

        # assign moltypes
        ligand.assign_moltype()
        complex.assign_moltype()

        return complex, ligand


# ==============
# Data/Libraries
# ==============
_dna_names = {
    'DA5_DT5': 'D5K',
    'DA5_DC5': 'D5L',
    'DA5_DG5': 'D5M',
    'DT5_DA5': 'D5N',
    'DT5_DC5': 'D5O',
    'DT5_DG5': 'D5P',
    'DC5_DA5': 'D5R',
    'DC5_DT5': 'D5S',
    'DC5_DG5': 'D5T',
    'DG5_DA5': 'D5X',
    'DG5_DT5': 'D5Y',
    'DG5_DC5': 'D5Z',
    'DA3_DT3': 'D3K',
    'DA3_DC3': 'D3L',
    'DA3_DG3': 'D3M',
    'DT3_DA3': 'D3N',
    'DT3_DC3': 'D3O',
    'DT3_DG3': 'D3P',
    'DC3_DA3': 'D3R',
    'DC3_DT3': 'D3S',
    'DC3_DG3': 'D3T',
    'DG3_DA3': 'D3X',
    'DG3_DT3': 'D3Y',
    'DG3_DC3': 'D3Z',
    # False names to avoid an error
    'DG3_DG3': 'FOO',
    'DC3_DC3': 'FOO',
    'DA3_DA3': 'FOO',
    'DT3_DT3': 'FOO',
    'DG5_DG5': 'FOO',
    'DC5_DC5': 'FOO',
    'DA5_DA5': 'FOO',
    'DT5_DT5': 'FOO',
    }

_rna_names = {
    'RA5_RU5': 'R5K',
    'RA5_RC5': 'R5L',
    'RA5_RG5': 'R5M',
    'RU5_RA5': 'R5N',
    'RU5_RC5': 'R5O',
    'RU5_RG5': 'R5P',
    'RC5_RA5': 'R5R',
    'RC5_RU5': 'R5S',
    'RC5_RG5': 'R5T',
    'RG5_RA5': 'R5X',
    'RG5_RU5': 'R5Y',
    'RG5_RC5': 'R5Z',
    'RA3_RU3': 'R3K',
    'RA3_RC3': 'R3L',
    'RA3_RG3': 'R3M',
    'RU3_RA3': 'R3N',
    'RU3_RC3': 'R3O',
    'RU3_RG3': 'R3P',
    'RC3_RA3': 'R3R',
    'RC3_RU3': 'R3S',
    'RC3_RG3': 'R3T',
    'RG3_RA3': 'R3X',
    'RG3_RU3': 'R3Y',
    'RG3_RC3': 'R3Z',
    # False names to avoid an error
    'RG3_RG3': 'FOO',
    'RC3_RC3': 'FOO',
    'RA3_RA3': 'FOO',
    'RU3_RU3': 'FOO',
    'RG5_RG5': 'FOO',
    'RC5_RC5': 'FOO',
    'RA5_RA5': 'FOO',
    'RU5_RU5': 'FOO',
    }
