"""Reference-anchored locus alignments built from minimap2 ``cs`` tags.

Every window (from the reference genome) was aligned pairwise to each genome by
minimap2. The ``cs`` string of that hit tells, base by base, what the genome shows at
every window position (match, substitution, deletion) and where the genome carries
extra bases (insertions). Stacking these per-genome strings over the window
coordinates yields a multiple alignment without running an MSA program. Insertions
of different genomes at the same window position are left-aligned against each
other (an approximation that is irrelevant for primer placement, because every
position touched by any insertion is masked anyway).
"""

from __future__ import annotations

import re
from typing import Dict, List, Mapping, Tuple

from .seqio import revcomp

_CS_OP = re.compile(r"(:[0-9]+|\*[a-zA-Z][a-zA-Z]|[=+\-][A-Za-z]+)")


def parse_cs(cs: str) -> List[Tuple[str, str]]:
    """Split a cs string into (op, payload): (':', '12'), ('*', 'ag'), ('+', 'acg'), ('-', 't')."""
    ops = []
    pos = 0
    for m in _CS_OP.finditer(cs):
        if m.start() != pos:
            raise ValueError(f"Cannot parse cs string near offset {pos}: {cs[pos:pos+20]!r}")
        tok = m.group(0)
        ops.append((tok[0], tok[1:]))
        pos = m.end()
    if pos != len(cs):
        raise ValueError(f"Trailing garbage in cs string: {cs[pos:]!r}")
    return ops


def flip_cs(ops: List[Tuple[str, str]]) -> List[Tuple[str, str]]:
    """Convert a minus-strand cs (rev-comp query vs forward target) to forward query orientation."""
    flipped = []
    for op, payload in reversed(ops):
        if op == ":":
            flipped.append((op, payload))
        elif op == "*":
            flipped.append((op, revcomp(payload[0]) + revcomp(payload[1])))
        else:  # '+', '-', '='
            flipped.append((op, revcomp(payload)))
    return flipped


def hit_states(
    window_seq: str, qstart: int, qend: int, strand: str, cs: str
) -> Tuple[List[str], Dict[int, str]]:
    """What one genome shows at each window position.

    Returns (states, insertions): ``states[i]`` is the genome base at window position
    ``i`` ('-' if the genome lacks it, 'N' outside the aligned part), and
    ``insertions[i]`` are genome bases that sit between window positions i-1 and i.
    """
    n = len(window_seq)
    states = ["N"] * n
    insertions: Dict[int, str] = {}
    ops = parse_cs(cs)
    if strand == "-":
        ops = flip_cs(ops)
    q = qstart
    for op, payload in ops:
        if op == ":":
            k = int(payload)
            states[q : q + k] = window_seq[q : q + k]
            q += k
        elif op == "=":
            k = len(payload)
            states[q : q + k] = payload.upper()
            q += k
        elif op == "*":
            states[q] = payload[0].upper()  # genome (PAF target) base
            q += 1
        elif op == "+":  # window has bases the genome lacks
            k = len(payload)
            states[q : q + k] = ["-"] * k
            q += k
        elif op == "-":  # genome has extra bases
            insertions[q] = insertions.get(q, "") + payload.upper()
    if q != qend:
        raise ValueError(f"cs string length ({q - qstart}) does not match qend-qstart ({qend - qstart})")
    return states, insertions


def build_alignment(
    window_seq: str,
    per_genome: Mapping[str, Tuple[int, int, str, str]],
) -> Dict[str, str]:
    """Stack per-genome (qstart, qend, strand, cs) into an alignment {genome_id: aligned_seq}.

    Column order: for each window position i, first the insertion block (width = the
    longest insertion any genome has before i), then position i itself.
    """
    n = len(window_seq)
    parsed = {gid: hit_states(window_seq, *hit) for gid, hit in per_genome.items()}
    ins_width = [0] * (n + 1)
    for _states, ins in parsed.values():
        for pos, seq in ins.items():
            ins_width[pos] = max(ins_width[pos], len(seq))
    out: Dict[str, str] = {}
    for gid, (states, ins) in parsed.items():
        parts: List[str] = []
        for i in range(n + 1):
            if ins_width[i]:
                s = ins.get(i, "")
                parts.append(s + "-" * (ins_width[i] - len(s)))
            if i < n:
                parts.append(states[i])
        out[gid] = "".join(parts)
    return out
