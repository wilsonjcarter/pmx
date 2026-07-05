"""Friendly shorthand for common post-translational modifications.

The generative ``mutate`` operation turns any WT->target morph into a hybrid, so
a PTM is just ``mutate(wt, target)``.  These codes spare the user from
memorising force-field residue names:

======================  ===================================  =================
modification            code                                 residue (WT)
======================  ===================================  =================
phosphorylation         ``p1`` monoanion (-1) / ``p2`` (-2)  SER / THR / TYR
methylation             ``me1`` / ``me2`` / ``me3`` (mono/   LYS
                        di/tri)
======================  ===================================  =================

Codes are case-insensitive; the **digit is required** (a bare ``p`` would
collide with proline's one-letter code).  :func:`resolve` returns ``None`` for
anything that is not a PTM code, so callers fall back to treating the token as a
literal residue name.

Target names follow the CHARMM36m convention (``SP1/SP2``, ``MLY/M3L`` ...).
Other force fields may name these residues differently; the caller should verify
the resolved target exists in the chosen force field and surface a clear error
otherwise.
"""

from __future__ import annotations

import re
from typing import Optional, Tuple

# (WT residue) -> {code: target residue name}, per modification family.
_PHOSPHO = {
    "SER": {"p1": "SP1", "p2": "SP2"},     # phospho-serine
    "THR": {"p1": "TP1", "p2": "TP2"},     # phospho-threonine
    "TYR": {"p1": "YP1", "p2": "YP2"},     # phospho-tyrosine
}
_METHYL = {
    "LYS": {"me1": "MLYS", "me2": "MLY", "me3": "M3L"},  # mono/di/tri-methyl-Lys
}

# net charge change (target - WT) for each code -- phosphorylation deprotonates
# the added phosphate (-1 or -2); methylation keeps the +1 ammonium.
_DELTA_Q = {"p1": -1, "p2": -2, "me1": 0, "me2": 0, "me3": 0}

# digit is mandatory: bare 'p' collides with proline, bare 'me' is ambiguous.
_CODE_RE = re.compile(r"^(?:p[12]|me[123])$")

# every PTM target residue name (the single source of truth, no scattered dicts)
TARGET_RESNAMES = frozenset(
    tgt for table in (_PHOSPHO, _METHYL)
    for codes in table.values() for tgt in codes.values())


def is_ptm_code(token: str) -> bool:
    """True if *token* is a recognised PTM shorthand (case-insensitive)."""
    return bool(_CODE_RE.match(token.strip().lower()))


def resolve(wt_resname: str, code: str) -> Optional[Tuple[str, int]]:
    """Resolve a PTM code against a wild-type residue.

    Returns ``(target_resname, delta_charge)``, or ``None`` if *code* is not a
    PTM code (so the caller treats it as a literal residue name).

    Raises
    ------
    ValueError
        If *code* is a valid PTM code but is not defined for *wt_resname*
        (e.g. ``p2`` on a lysine).
    """
    c = code.strip().lower()
    if not _CODE_RE.match(c):
        return None
    wt = wt_resname.strip().upper()
    table = _PHOSPHO if c[0] == "p" else _METHYL
    family = "phosphorylation" if c[0] == "p" else "methylation"
    if wt not in table or c not in table[wt]:
        supported = "; ".join(
            f"{r} -> {'/'.join(k for k in tbl)}"
            for tbl_family in (_PHOSPHO, _METHYL)
            for r, tbl in tbl_family.items())
        raise ValueError(
            f"PTM code '{code}' ({family}) is not defined for residue "
            f"'{wt_resname}'.  Supported: {supported}.")
    return table[wt][c], _DELTA_Q[c]


def target_for(wt_resname: str, token: str) -> Optional[str]:
    """Target residue for *token* (a PTM code *or* a direct PTM target name),
    or ``None`` if *token* is neither -- so callers route PTMs through the
    generative engine and fall back to their normal path otherwise.

    This is the single dispatch point that lets ``mutate`` recognise a
    phosphorylation/methylation and hand it to the generative builder instead of
    a hand-authored library entry.
    """
    if is_ptm_code(token):
        r = resolve(wt_resname, token)          # may raise for an ill-matched WT
        return r[0] if r else None
    if token.strip().upper() in TARGET_RESNAMES:
        return token.strip().upper()
    return None
