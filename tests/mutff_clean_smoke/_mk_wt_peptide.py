#!/usr/bin/env python
"""Build a clean single-state WT capped peptide for the PTM smoke matrix.

Usage:
    _mk_wt_peptide.py <template.pdb> <WT-resname> <out.pdb> <ff-name> [resid]

The template is ``ACE-ALA-CYM-ALA-NME`` (resid 3 = CYM).  We replace resid 3 with
a clean WT residue (SER/THR/TYR/LYS/...) built *template-free* from the force
field RTP (``TargetResidue.from_rtp``) and grafted onto the CYM backbone frame
with the very same Kabsch superposition the generative engine uses
(``_superpose_target``).  The result is a fully protonated, single-state WT
peptide -- a valid starting point for the subsequent generative WT->PTM mutation
(the smoke test's actual subject).

Why not a plain CYM->WT ``pmx mutate``?  That emits a dual-state *hybrid*; here we
want a pristine WT end-state, so we graft the built residue in directly.
"""
import sys

from pmx.model import Model
from pmx.covalent.mutate import (_load_rtp, _assign_rtp_entries,
                                 _build_target_residue, _superpose_target)
from pmx.covalent.mapping import map_residues, attach_bonds_from_rtp
from pmx.covalent.mapping import _element_of  # noqa: F401  (kept for parity)


def main():
    if len(sys.argv) < 5:
        sys.exit(__doc__)
    template, wt, out, ff = sys.argv[1:5]
    wt = wt.upper()
    resid = int(sys.argv[5]) if len(sys.argv) > 5 else 3

    m = Model(template, bPDBTER=True, rename_atoms=True, scale_coords='A')
    if m.unity == "nm":
        m.nm2a()

    rtp, _ = _load_rtp(ff)
    residue1 = m.fetch_residue(idx=resid)          # CYM at resid 3
    chain = residue1.chain
    term = ("C" if residue1 == chain.cterminus()
            else ("N" if residue1 == chain.nterminus() else None))

    # parameterize + bond the WT (CYM) so the graph mapping has atom types/bonds
    _assign_rtp_entries(residue1, rtp, term, ff=ff)
    if not any(getattr(a, "bonds", None) for a in residue1.atoms):
        attach_bonds_from_rtp(residue1, ff)

    # build the clean target WT residue from the RTP (params + embedded geometry)
    residue2 = _build_target_residue(wt, ff, rtp, term, residue1)
    residue2.set_resname(wt)

    # graft: place the built residue in the CYM backbone frame.  Fit on the
    # backbone atoms ONLY (not the whole residue) so the flanking peptide bonds
    # stay at their true length -- a whole-residue Kabsch compromises the N/C
    # position and stretches the attachment bond.
    mapping = map_residues(residue1, residue2)
    pairs = list(zip(mapping.pairs1(), mapping.pairs2()))
    _BB = {"N", "CA", "C", "O"}
    bb_pairs = [(i, j) for i, j in pairs if residue1.atoms[i].name in _BB]
    _superpose_target(residue1, residue2, bb_pairs if len(bb_pairs) >= 3 else pairs)

    # swap the clean residue into the chain at resid 3, then drop the CYM
    pos = chain.residues.index(residue1)
    chain.insert_residue(pos, residue2, newResNum=resid)
    m.remove_residue(residue1)

    m.write(out)
    print(f"log_> _mk_wt_peptide: {template} resid {resid} CYM -> {wt} "
          f"({len(residue2.atoms)} atoms) -> {out}")


if __name__ == "__main__":
    main()
