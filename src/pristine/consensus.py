"""Consensus and primer-template construction from a locus alignment."""

from __future__ import annotations

from typing import Dict, Iterable, List, Sequence

IUPAC_CODES = {
    frozenset("A"): "A", frozenset("C"): "C", frozenset("G"): "G", frozenset("T"): "T",
    frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("GC"): "S", frozenset("AT"): "W",
    frozenset("GT"): "K", frozenset("AC"): "M", frozenset("ACG"): "V", frozenset("ACT"): "H",
    frozenset("AGT"): "D", frozenset("CGT"): "B", frozenset("ACGT"): "N",
}


def _check(seqs: Sequence[str]) -> int:
    if not seqs:
        raise ValueError("No sequences")
    n = len(seqs[0])
    if any(len(s) != n for s in seqs):
        raise ValueError("Aligned sequences differ in length")
    return n


def iupac_consensus(seqs: Sequence[str]) -> str:
    """IUPAC consensus ignoring gaps (documentation output, not used for design)."""
    n = _check(seqs)
    out: List[str] = []
    for i in range(n):
        bases = {s[i].upper() for s in seqs if s[i] not in "-" and s[i].upper() in "ACGT"}
        out.append(IUPAC_CODES.get(frozenset(bases), "N") if bases else "N")
    return "".join(out)


def masked_template(seqs: Sequence[str]) -> str:
    """Primer design template: a base where ALL sequences agree (A/C/G/T, no gap), else N.

    Length equals the alignment length so alignment coordinates can be used directly
    for SEQUENCE_TARGET and for locating informative positions in amplicons.
    """
    n = _check(seqs)
    out: List[str] = []
    for i in range(n):
        col = {s[i].upper() for s in seqs}
        out.append(col.pop() if len(col) == 1 and col.issubset(set("ACGT")) and col else "N")
    return "".join(out)
