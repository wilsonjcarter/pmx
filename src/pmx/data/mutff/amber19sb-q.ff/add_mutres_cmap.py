#!/usr/bin/env python3
"""
Add CMAP entries to amber19sb-q.ff for pmx hybrid residues in mutres.rtp.

For ff19SB-style residue-specific CMAPs, pdb2gmx (>=2026) constructs the
CMAP atom-type lookup as "<type>-<residue_name>" (e.g. "XC-ALA" for ALA's
CA). Hybrid residues like "A2R" (Ala->Arg) have residue name A2R, so the
lookup looks for "XC-A2R", which doesn't exist.

For each hybrid residue:
  1. Append a [ cmap ] section to the mutres.rtp entry, mirroring the
     aminoacids.rtp form: `-C N CA C +N`.
  2. Append a new cmaptypes block to cmap.itp keyed by the hybrid name
     (e.g. "C-* N-A2R XC-A2R C-A2R N-*") containing the A-state residue's
     CMAP grid (e.g. ALA's grid for A2R).

The A-state is the first letter of the hybrid name (e.g. "A2R" -> A -> ALA).
"""
import re, os, shutil, sys

ROOT = '/sessions/sharp-determined-shannon/mnt/extra_ff/amber19sb-q.ff'
RTP  = os.path.join(ROOT, 'mutres.rtp')
CMAP = os.path.join(ROOT, 'cmap.itp')

# Single-letter A-state code -> full residue name
A_STATE_MAP = {
    'A':'ALA','R':'ARG','B':'ASH','N':'ASN','D':'ASP','C':'CYS',
    'J':'GLH','Q':'GLN','E':'GLU','G':'GLY','H':'HID','X':'HIE','Z':'HIP',
    'I':'ILE','L':'LEU','O':'LYN','K':'LYS','M':'MET','F':'PHE','P':'PRO',
    'S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL',
}

# ----------------------------------------------------------------------
# Step 1: Parse existing cmap.itp -> dict of residue_name -> full block text.
# Each block starts with "C-* N-RES XC-RES C-RES N-* 1 24 24\" line and runs
# until the next "C-*" line at column 0 (or EOF).
# ----------------------------------------------------------------------
def parse_cmap_blocks(cmap_text):
    """Return list of (residue_name, central_ca_type, raw_block_text) for each entry,
       and the index where the [ cmaptypes ] section's last block ends so we know
       where to append new blocks.
    """
    blocks = []
    # Identify start of each cmap entry
    pattern = re.compile(r'^(C-\*\s+N-(\w+)\s+(XC|CA)-\w+\s+C-\w+\s+N-\*[^\n]*)$', re.MULTILINE)
    matches = list(pattern.finditer(cmap_text))
    for i, m in enumerate(matches):
        start = m.start()
        end = matches[i+1].start() if i+1 < len(matches) else len(cmap_text)
        block = cmap_text[start:end]
        blocks.append({'name': m.group(2), 'ca_type': m.group(3), 'text': block.rstrip()+'\n'})
    return blocks

with open(CMAP) as f: cmap_txt = f.read()
blocks = parse_cmap_blocks(cmap_txt)
res_to_block = {b['name']: b for b in blocks}
print(f'Parsed {len(blocks)} cmap blocks: {sorted(res_to_block)}')

# ----------------------------------------------------------------------
# Step 2: Get list of hybrid residue names from mutres.rtp
# ----------------------------------------------------------------------
with open(RTP) as f: rtp_txt = f.read()
hybrids = re.findall(r'^\[ ([A-Z]2[A-Z]) \]', rtp_txt, re.MULTILINE)
print(f'Hybrid residues to process: {len(hybrids)}')

# ----------------------------------------------------------------------
# Step 3: For each hybrid, build a new cmaptypes block from the A-state grid.
# ----------------------------------------------------------------------
new_cmap_blocks = []
unhandled = []
for hybrid in hybrids:
    a_letter = hybrid[0]
    a_res = A_STATE_MAP.get(a_letter)
    if a_res is None or a_res not in res_to_block:
        unhandled.append((hybrid, a_letter, a_res))
        continue
    base_block = res_to_block[a_res]
    # Replace the residue suffix in the header line
    base_text = base_block['text']
    ca_type = base_block['ca_type']  # XC for most, CA for CNX
    header_re = re.compile(rf'^C-\*\s+N-{a_res}\s+{ca_type}-{a_res}\s+C-{a_res}\s+N-\*([^\n]*)$', re.MULTILINE)
    new_header_match = header_re.search(base_text)
    if not new_header_match:
        print(f'  WARN: header lookup failed for {hybrid} (A-state {a_res})')
        unhandled.append((hybrid, a_letter, a_res))
        continue
    new_header = f'C-* N-{hybrid} {ca_type}-{hybrid} C-{hybrid} N-*' + new_header_match.group(1)
    new_block_text = base_text[:new_header_match.start()] + new_header + base_text[new_header_match.end():]
    # Prefix with a comment line for traceability
    new_cmap_blocks.append(f'; --- pmx hybrid {hybrid} (A-state {a_res}) ---\n{new_block_text}')

print(f'Built {len(new_cmap_blocks)} new cmap blocks; {len(unhandled)} unhandled')
if unhandled[:5]: print('  unhandled examples:', unhandled[:5])

# ----------------------------------------------------------------------
# Step 4: Append new cmap blocks to cmap.itp
# ----------------------------------------------------------------------
new_cmap_txt = cmap_txt.rstrip() + '\n\n' + \
    '; ====================================================================\n' + \
    '; pmx HYBRID RESIDUE CMAP entries (auto-generated)\n' + \
    '; Each hybrid residue (e.g. A2R = Ala->Arg) inherits the A-state grid.\n' + \
    '; ====================================================================\n\n' + \
    ''.join(new_cmap_blocks)
with open(CMAP, 'w') as f: f.write(new_cmap_txt)
print(f'Wrote {CMAP} (added {len(new_cmap_blocks)} blocks)')

# ----------------------------------------------------------------------
# Step 5: Add [ cmap ] section to each mutres.rtp entry.
# Insert it after the LAST [ subsection ] of each hybrid residue block.
# ----------------------------------------------------------------------
def add_cmap_to_mutres(rtp_text, hybrid):
    """Append a [ cmap ] section to the end of [ hybrid ] residue block."""
    m = re.search(rf'^\[ {hybrid} \]', rtp_text, re.MULTILINE)
    if not m: return rtp_text, False
    # Find end of this residue block: next "[ XYZ ]" header at line start, or EOF
    after = rtp_text[m.end():]
    nxt = re.search(r'^\[ [A-Z][A-Z0-9]+ \]', after, re.MULTILINE)
    end_of_block = m.end() + (nxt.start() if nxt else len(after))
    block = rtp_text[m.end():end_of_block]
    # Skip if [ cmap ] already present
    if re.search(r'^\s*\[ cmap \]', block, re.MULTILINE):
        return rtp_text, False
    # Build the cmap section text
    cmap_section = '\n [ cmap ]\n -C    N     CA    C     +N\n'
    # Insert just before end_of_block, preserving any trailing blank lines
    # The block usually ends with "[ dihedrals ]" possibly empty + blank lines
    # before next header. We insert cmap_section right before the next header.
    new_block = block.rstrip() + cmap_section + '\n'
    new_rtp = rtp_text[:m.end()] + new_block + rtp_text[end_of_block:]
    return new_rtp, True

with open(RTP) as f: rtp_txt = f.read()
modified = 0
for hybrid in hybrids:
    rtp_txt, did = add_cmap_to_mutres(rtp_txt, hybrid)
    if did: modified += 1
with open(RTP,'w') as f: f.write(rtp_txt)
print(f'Added [ cmap ] to {modified} hybrid residues in mutres.rtp')

print('DONE')
