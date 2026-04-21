#!/usr/bin/env python
"""
Generate mutres_term.rtp and mutres_term.mtp files for terminal residue deletions.

For C-terminal deletion (XdeC):
  - The n-1 residue gains a dummy OT2 atom (DUM_OC, charge 0.0 in state A)
  - In state B, the residue becomes a true C-terminus:
    * C morphs from C -> CC, charge 0.51 -> 0.34
    * O morphs from O -> OC, charge -0.51 -> -0.67
    * DOT2 morphs from DUM_OC -> OC, charge 0.0 -> -0.67
  - The +N bond and related impropers are removed (no next residue in state B)

For N-terminal deletion (XdeN):
  - The n+1 residue gains two dummy protons (DH2, DH3; DUM_HC, charge 0.0 in state A)
  - In state B, the residue becomes a true N-terminus:
    * N morphs from NH1 -> NH3, charge -0.47 -> -0.30
    * HN morphs from H -> HC, name stays HN in topology but type changes, charge 0.31 -> 0.33
    * DH2 morphs from DUM_HC -> HC, charge 0.0 -> 0.33
    * DH3 morphs from DUM_HC -> HC, charge 0.0 -> 0.33
    * CA type may change (CT1 -> CT1 stays, but charge 0.07 -> 0.21)
    * HA type may change (HB1 -> HB1, charge 0.09 -> 0.10)
  - The -C bond and related impropers are removed (no previous residue in state B)

Special cases:
  - GLY: CA is CT2 (not CT1), HA1/HA2 (not HA), different N-term charges
  - PRO: No HN atom, N-terminal patch is PRO-NH2+ (NP type, different charges)
         For N-term deletion neighbor, PRO gains 1 dummy proton (not 2)
"""

import os
import re


FF_DIR = None  # Set at runtime

THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
    'CYS': 'C', 'CYM': 'CM', 'GLN': 'Q', 'GLU': 'E',
    'GLY': 'G', 'HSD': 'X', 'HSE': 'H', 'HSP': 'Z',
    'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M',
    'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T',
    'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}

NTER_STANDARD = {
    'N':  {'type': 'NH3', 'charge': -0.30},
    'CA': {'type': 'CT1', 'charge':  0.21},
    'HA': {'type': 'HB1', 'charge':  0.10},
    'HN': 'delete',
    'add': [('H1', 'HC', 0.33), ('H2', 'HC', 0.33), ('H3', 'HC', 0.33)],
}

NTER_GLY = {
    'N':   {'type': 'NH3', 'charge': -0.30},
    'CA':  {'type': 'CT2', 'charge':  0.13},
    'HN': 'delete',
    'add': [('H1', 'HC', 0.33), ('H2', 'HC', 0.33), ('H3', 'HC', 0.33)],
}

NTER_PRO = {
    'N':  {'type': 'NP',  'charge': -0.07},
    'CA': {'type': 'CP1', 'charge':  0.16},
    'CD': {'type': 'CP3', 'charge':  0.16},
    'add': [('HN1', 'HC', 0.24), ('HN2', 'HC', 0.24)],
}

CTER = {
    'C':  {'type': 'CC', 'charge':  0.34},
    'O':  {'type': 'OC', 'charge': -0.67, 'rename': 'OT1'},
    'add': [('OT2', 'OC', -0.67)],
}


def parse_residue_from_rtp(rtp_path, resname):
    """Parse a single residue entry from an .rtp file."""
    atoms = []
    bonds = []
    impropers = []
    cmap = []

    with open(rtp_path) as f:
        lines = f.readlines()

    in_residue = False
    current_section = None

    for line in lines:
        stripped = line.strip()

        if re.match(r'^\[\s*' + re.escape(resname) + r'\s*\]', stripped):
            in_residue = True
            current_section = None
            continue

        if in_residue:
            if re.match(r'^\[\s*[A-Z]', stripped) and not re.match(r'^\[\s*(atoms|bonds|impropers|dihedrals|cmap)\s*\]', stripped, re.IGNORECASE):
                break

            m = re.match(r'^\[\s*(atoms|bonds|impropers|dihedrals|cmap)\s*\]', stripped, re.IGNORECASE)
            if m:
                current_section = m.group(1).lower()
                continue

            if not stripped or stripped.startswith(';'):
                continue

            if current_section == 'atoms':
                parts = stripped.split()
                if len(parts) >= 4:
                    atoms.append({
                        'name': parts[0],
                        'type': parts[1],
                        'charge': float(parts[2]),
                        'index': int(parts[3]),
                    })

            elif current_section == 'bonds':
                parts = stripped.split(';')[0].split()
                if len(parts) >= 2:
                    bonds.append((parts[0], parts[1]))

            elif current_section == 'impropers':
                parts = stripped.split(';')[0].split()
                if len(parts) >= 4:
                    impropers.append((parts[0], parts[1], parts[2], parts[3]))

            elif current_section == 'cmap':
                parts = stripped.split(';')[0].split()
                if len(parts) >= 5:
                    cmap.append(tuple(parts[:5]))

    return {
        'atoms': atoms,
        'bonds': bonds,
        'impropers': impropers,
        'cmap': cmap,
    }


def generate_cterm_neighbor_rtp(resname, res_data):
    """Generate .rtp entry for n-1 residue (C-terminal deletion neighbor)."""
    one = THREE_TO_ONE[resname]
    hybrid_name = f"{one}deC"

    atoms = list(res_data['atoms'])
    bonds = list(res_data['bonds'])
    impropers = list(res_data['impropers'])
    cmap = list(res_data['cmap'])

    next_idx = max(a['index'] for a in atoms) + 1
    atoms.append({
        'name': 'DOT2',
        'type': 'DUM_OC',
        'charge': 0.00000,
        'index': next_idx,
    })

    bonds.append(('DOT2', 'C'))
    impropers.append(('C', 'CA', 'DOT2', 'O'))

    lines = []
    lines.append(f"[ {hybrid_name} ] ; {resname} as C-term deletion neighbor\n")
    lines.append("")
    lines.append(" [ atoms ]")
    for a in atoms:
        lines.append(f"  {a['name']:>5s}   {a['type']:<16s} {a['charge']:8.5f}  {a['index']}")
    lines.append("")
    lines.append(" [ bonds ]")
    for b in bonds:
        lines.append(f"  {b[0]:>5s}  {b[1]:>5s}")
    lines.append("")
    lines.append(" [ impropers ]")
    for imp in impropers:
        lines.append(f"  {imp[0]:>5s}  {imp[1]:>5s}  {imp[2]:>5s}  {imp[3]:>5s}")
    if cmap:
        lines.append("")
        lines.append(" [ cmap ]")
        for cm in cmap:
            lines.append(f"  {cm[0]:>5s}  {cm[1]:>5s}  {cm[2]:>5s}  {cm[3]:>5s}  {cm[4]:>5s}")
    lines.append("")
    return "\n".join(lines), hybrid_name


def generate_nterm_neighbor_rtp(resname, res_data):
    """Generate .rtp entry for n+1 residue (N-terminal deletion neighbor)."""
    one = THREE_TO_ONE[resname]
    hybrid_name = f"{one}deN"

    atoms = list(res_data['atoms'])
    bonds = list(res_data['bonds'])
    impropers = list(res_data['impropers'])
    cmap = list(res_data['cmap'])

    next_idx = max(a['index'] for a in atoms) + 1

    if resname == 'PRO':
        atoms.append({'name': 'DHN1', 'type': 'DUM_HC', 'charge': 0.00000, 'index': next_idx})
        atoms.append({'name': 'DHN2', 'type': 'DUM_HC', 'charge': 0.00000, 'index': next_idx + 1})
        bonds.append(('DHN1', 'N'))
        bonds.append(('DHN2', 'N'))
    else:
        atoms.append({'name': 'DH2', 'type': 'DUM_HC', 'charge': 0.00000, 'index': next_idx})
        atoms.append({'name': 'DH3', 'type': 'DUM_HC', 'charge': 0.00000, 'index': next_idx + 1})
        bonds.append(('DH2', 'N'))
        bonds.append(('DH3', 'N'))

    lines = []
    lines.append(f"[ {hybrid_name} ] ; {resname} as N-term deletion neighbor\n")
    lines.append("")
    lines.append(" [ atoms ]")
    for a in atoms:
        lines.append(f"  {a['name']:>5s}   {a['type']:<16s} {a['charge']:8.5f}  {a['index']}")
    lines.append("")
    lines.append(" [ bonds ]")
    for b in bonds:
        lines.append(f"  {b[0]:>5s}  {b[1]:>5s}")
    lines.append("")
    lines.append(" [ impropers ]")
    for imp in impropers:
        lines.append(f"  {imp[0]:>5s}  {imp[1]:>5s}  {imp[2]:>5s}  {imp[3]:>5s}")
    if cmap:
        lines.append("")
        lines.append(" [ cmap ]")
        for cm in cmap:
            lines.append(f"  {cm[0]:>5s}  {cm[1]:>5s}  {cm[2]:>5s}  {cm[3]:>5s}  {cm[4]:>5s}")
    lines.append("")
    return "\n".join(lines), hybrid_name


def generate_cterm_neighbor_mtp(resname, res_data):
    """Generate .mtp entry for n-1 residue (C-terminal deletion neighbor).

    State A: normal mid-chain residue
    State B: C-terminus (C->CC, O->OC renamed OT1, DOT2->OC as OT2)
    """
    one = THREE_TO_ONE[resname]
    hybrid_name = f"{one}deC"
    atoms_a = res_data['atoms']

    lines = []
    lines.append(f"[ {hybrid_name} ] ; {resname} as C-term deletion neighbor\n")
    lines.append("")

    lines.append(" [ morphes ]")
    for a in atoms_a:
        name = a['name']
        typeA = a['type']
        if name == 'C':
            lines.append(f"  {name:>5s}  {typeA:>10s} ->  {name:>5s}  {'CC':>10s}")
        elif name == 'O':
            lines.append(f"  {name:>5s}  {typeA:>10s} ->  {'OT1':>5s}  {'OC':>10s}")
        else:
            lines.append(f"  {name:>5s}  {typeA:>10s} ->  {name:>5s}  {typeA:>10s}")
    lines.append(f"  {'DOT2':>5s}  {'DUM_OC':>10s} ->  {'OT2':>5s}  {'OC':>10s}")
    lines.append("")

    lines.append(" [ atoms ]")
    for a in atoms_a:
        name = a['name']
        typeA = a['type']
        chargeA = a['charge']
        massA = _get_mass(typeA)
        if name == 'C':
            typeB, chargeB, massB = 'CC', 0.34, 12.011
        elif name == 'O':
            typeB, chargeB, massB = 'OC', -0.67, 15.999
        else:
            typeB, chargeB, massB = typeA, chargeA, massA
        type_eq = 'types ==' if typeA == typeB else 'types !='
        charge_eq = 'charge ==' if abs(chargeA - chargeB) < 1e-6 else 'charge !='
        lines.append(f"  {name:>5s}  {typeA:>10s}  {chargeA:10.6f}      1  {massA:10.6f}  {typeB:>10s}  {chargeB:10.6f}  {massB:10.6f}   ;  {type_eq} | {charge_eq}")
    lines.append(f"  {'DOT2':>5s}  {'DUM_OC':>10s}  {0.0:10.6f}      1  {15.999:10.6f}  {'OC':>10s}  {-0.67:10.6f}  {15.999:10.6f}   ;  types != | charge !=")
    lines.append("")

    lines.append(" [ coords ]")
    for _ in atoms_a:
        lines.append("    0.000    0.000    0.000")
    lines.append("    0.000    0.000    0.000")  # DOT2
    lines.append("")

    lines.append(" [ impropers ]")
    lines.append(f"  {'N':>5s}  {'-C':>5s}  {'CA':>5s}  {'HN':>5s}     default-A                 default-B")
    lines.append(f"  {'C':>5s}  {'CA':>5s}  {'DOT2':>5s}  {'O':>5s}     default-A                 default-B")
    lines.append("")
    lines.append(" [ dihedrals ]")
    lines.append("")
    lines.append(" [ rotations ]")
    lines.append("")

    return "\n".join(lines), hybrid_name


def generate_nterm_neighbor_mtp(resname, res_data):
    """Generate .mtp entry for n+1 residue (N-terminal deletion neighbor).

    State A: normal mid-chain residue
    State B: N-terminus (N->NH3, HN->HC as H1, DH2->HC, DH3->HC)
    """
    one = THREE_TO_ONE[resname]
    hybrid_name = f"{one}deN"
    atoms_a = res_data['atoms']

    if resname == 'PRO':
        nter = NTER_PRO
    elif resname == 'GLY':
        nter = NTER_GLY
    else:
        nter = NTER_STANDARD

    lines = []
    lines.append(f"[ {hybrid_name} ] ; {resname} as N-term deletion neighbor\n")
    lines.append("")

    lines.append(" [ morphes ]")
    for a in atoms_a:
        name = a['name']
        typeA = a['type']
        if name == 'N':
            typeB = nter['N']['type']
            lines.append(f"  {name:>5s}  {typeA:>10s} ->  {name:>5s}  {typeB:>10s}")
        elif name == 'HN' and resname != 'PRO':
            lines.append(f"  {'HN':>5s}  {typeA:>10s} ->  {'H1':>5s}  {'HC':>10s}")
        elif name == 'CA':
            typeB = nter.get('CA', {}).get('type', typeA)
            lines.append(f"  {name:>5s}  {typeA:>10s} ->  {name:>5s}  {typeB:>10s}")
        elif name == 'HA' and 'HA' in nter:
            typeB = nter['HA']['type']
            lines.append(f"  {name:>5s}  {typeA:>10s} ->  {name:>5s}  {typeB:>10s}")
        elif name == 'CD' and resname == 'PRO' and 'CD' in nter:
            typeB = nter['CD']['type']
            lines.append(f"  {name:>5s}  {typeA:>10s} ->  {name:>5s}  {typeB:>10s}")
        else:
            lines.append(f"  {name:>5s}  {typeA:>10s} ->  {name:>5s}  {typeA:>10s}")
    if resname == 'PRO':
        lines.append(f"  {'DHN1':>5s}  {'DUM_HC':>10s} ->  {'HN1':>5s}  {'HC':>10s}")
        lines.append(f"  {'DHN2':>5s}  {'DUM_HC':>10s} ->  {'HN2':>5s}  {'HC':>10s}")
    else:
        lines.append(f"  {'DH2':>5s}  {'DUM_HC':>10s} ->  {'H2':>5s}  {'HC':>10s}")
        lines.append(f"  {'DH3':>5s}  {'DUM_HC':>10s} ->  {'H3':>5s}  {'HC':>10s}")
    lines.append("")

    lines.append(" [ atoms ]")
    for a in atoms_a:
        name = a['name']
        typeA = a['type']
        chargeA = a['charge']
        massA = _get_mass(typeA)
        if name in nter and nter[name] != 'delete':
            typeB = nter[name]['type']
            chargeB = nter[name]['charge']
            massB = _get_mass(typeB)
        elif name == 'HN' and resname != 'PRO':
            typeB, chargeB, massB = 'HC', 0.33, 1.008
        else:
            typeB, chargeB, massB = typeA, chargeA, massA
        type_eq = 'types ==' if typeA == typeB else 'types !='
        charge_eq = 'charge ==' if abs(chargeA - chargeB) < 1e-6 else 'charge !='
        lines.append(f"  {name:>5s}  {typeA:>10s}  {chargeA:10.6f}      1  {massA:10.6f}  {typeB:>10s}  {chargeB:10.6f}  {massB:10.6f}   ;  {type_eq} | {charge_eq}")
    if resname == 'PRO':
        for dname in ['DHN1', 'DHN2']:
            lines.append(f"  {dname:>5s}  {'DUM_HC':>10s}  {0.0:10.6f}      1  {1.008:10.6f}  {'HC':>10s}  {0.24:10.6f}  {1.008:10.6f}   ;  types != | charge !=")
    else:
        for dname in ['DH2', 'DH3']:
            lines.append(f"  {dname:>5s}  {'DUM_HC':>10s}  {0.0:10.6f}      1  {1.008:10.6f}  {'HC':>10s}  {0.33:10.6f}  {1.008:10.6f}   ;  types != | charge !=")
    lines.append("")

    lines.append(" [ coords ]")
    for _ in atoms_a:
        lines.append("    0.000    0.000    0.000")
    if resname == 'PRO':
        lines.append("    0.000    0.000    0.000")  # DHN1
        lines.append("    0.000    0.000    0.000")  # DHN2
    else:
        lines.append("    0.000    0.000    0.000")  # DH2
        lines.append("    0.000    0.000    0.000")  # DH3
    lines.append("")

    lines.append(" [ impropers ]")
    if resname == 'PRO':
        lines.append(f"  {'N':>5s}  {'-C':>5s}  {'CA':>5s}  {'CD':>5s}     default-A                 default-B")
    else:
        lines.append(f"  {'N':>5s}  {'-C':>5s}  {'CA':>5s}  {'HN':>5s}     default-A                 default-B")
    lines.append(f"  {'C':>5s}  {'CA':>5s}  {'+N':>5s}  {'O':>5s}     default-A                 default-B")
    lines.append("")
    lines.append(" [ dihedrals ]")
    lines.append("")
    lines.append(" [ rotations ]")
    lines.append("")

    return "\n".join(lines), hybrid_name


def _get_mass(atom_type):
    """Return mass based on CHARMM atom type."""
    type_upper = atom_type.upper().replace('DUM_', '')
    if type_upper.startswith('H') or type_upper in ('HC', 'H', 'HS', 'HB1', 'HB2', 'HA1', 'HA2', 'HA3', 'HP', 'HGA3'):
        return 1.008
    elif type_upper.startswith('C') or type_upper in ('CC', 'CD', 'CT1', 'CT2', 'CT3', 'CP1', 'CP2', 'CP3', 'CA'):
        return 12.011
    elif type_upper.startswith('N') or type_upper in ('NH1', 'NH2', 'NH3', 'NP', 'NC2'):
        return 14.007
    elif type_upper.startswith('O') or type_upper in ('OC', 'OH1', 'OB', 'O', 'OS'):
        return 15.999
    elif type_upper.startswith('S') or type_upper in ('S', 'SM'):
        return 32.060
    else:
        print(f"  WARNING: Unknown mass for type {atom_type}, defaulting to 12.011")
        return 12.011


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Generate mutres_term.rtp and mutres_term.mtp for terminal residue deletions")
    parser.add_argument('--ffdir', required=True,
                        help='Path to the force field directory (e.g., charmm36m-mut.ff)')
    parser.add_argument('--outdir', default=None,
                        help='Output directory (default: same as ffdir)')
    args = parser.parse_args()

    ffdir = args.ffdir
    outdir = args.outdir or ffdir
    rtp_path = os.path.join(ffdir, 'merged.rtp')

    if not os.path.exists(rtp_path):
        print(f"ERROR: {rtp_path} not found")
        return

    residues = [
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'CYM', 'GLN', 'GLU',
        'GLY', 'HSD', 'HSE', 'HSP', 'ILE', 'LEU', 'LYS', 'MET',
        'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    ]

    res_data = {}
    for resname in residues:
        data = parse_residue_from_rtp(rtp_path, resname)
        if not data['atoms']:
            print(f"WARNING: Could not parse {resname} from {rtp_path}")
            continue
        res_data[resname] = data

    rtp_lines = []
    rtp_lines.append("[ bondedtypes ]")
    rtp_lines.append("; Col 1: Type of bond")
    rtp_lines.append("; Col 2: Type of angles")
    rtp_lines.append("; Col 3: Type of proper dihedrals")
    rtp_lines.append("; Col 4: Type of improper dihedrals")
    rtp_lines.append("; Col 5: Generate all dihedrals if 1, only heavy atoms of 0.")
    rtp_lines.append("; Col 6: Number of excluded neighbors for nonbonded interactions")
    rtp_lines.append("; Col 7: Generate 1,4 interactions between pairs of hydrogens if 1")
    rtp_lines.append("; Col 8: Remove impropers over the same bond as a proper if it is 1")
    rtp_lines.append("; bonds  angles  dihedrals  impropers all_dihedrals nrexcl HH14 RemoveDih")
    rtp_lines.append("     1       5          9          2        1         3      1     0")
    rtp_lines.append("")

    mtp_lines = []

    for resname in sorted(res_data.keys()):
        data = res_data[resname]

        rtp_entry, name_c = generate_cterm_neighbor_rtp(resname, data)
        rtp_lines.append(rtp_entry)
        mtp_entry, _ = generate_cterm_neighbor_mtp(resname, data)
        mtp_lines.append(mtp_entry)

        rtp_entry, name_n = generate_nterm_neighbor_rtp(resname, data)
        rtp_lines.append(rtp_entry)
        mtp_entry, _ = generate_nterm_neighbor_mtp(resname, data)
        mtp_lines.append(mtp_entry)

        print(f"  Generated: {name_c}, {name_n}")

    rtp_out = os.path.join(outdir, 'mutres_term.rtp')
    mtp_out = os.path.join(outdir, 'mutres_term.mtp')

    with open(rtp_out, 'w') as f:
        f.write("\n".join(rtp_lines))
    print(f"\nWrote {rtp_out}")

    with open(mtp_out, 'w') as f:
        f.write("\n".join(mtp_lines))
    print(f"Wrote {mtp_out}")


if __name__ == '__main__':
    main()
