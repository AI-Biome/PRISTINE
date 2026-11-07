"""Cross-pair check of the selected primer pairs (for pooling them in one reaction).

Each pair is screened on its own by ``pipeline.evaluate_pair`` (F/R, F/F and R/R of the
same pair). When several pairs are pooled, a primer of one pair can also prime together
with a primer of another pair, for example when two selected loci lie close together
on the genome. This module enumerates those cross-pair products in every genome with
the same binding model as the in silico PCR.
"""

from __future__ import annotations

from bisect import bisect_right
from typing import List, Mapping, Sequence, Tuple

import pandas as pd

from .insilico_pcr import BindingParams, GenomeIndex, find_sites

COLUMNS = ["Genome", "Contig", "Start", "End", "Size", "Forward_Primer", "Reverse_Primer",
           "Forward_Locus", "Reverse_Locus", "Forward_Tm", "Reverse_Tm"]


def _primer_names(pairs: Sequence[Tuple[str, str, str]]) -> List[Tuple[str, str, str]]:
    """(name, locus, sequence) for every primer of every (locus, left, right) pair."""
    out = []
    for locus, left, right in pairs:
        out.append((f"{locus}_LEFT", locus, left))
        out.append((f"{locus}_RIGHT", locus, right))
    return out


def cross_pair_products(
    pairs: Sequence[Tuple[str, str, str]],
    indexes: Mapping[str, GenomeIndex],
    binding: BindingParams,
    max_product_size: int,
) -> pd.DataFrame:
    """Products formed by primers of two different pairs, in every genome.

    `pairs` are (locus, left_primer, right_primer). A product is a plus-strand site and
    a downstream minus-strand site at most `max_product_size` apart whose primers belong
    to different loci. Products within one pair are not reported here; the per-pair in
    silico PCR already covers them.
    """
    primers = _primer_names(pairs)
    rows = []
    if len(pairs) < 2:
        return pd.DataFrame(rows, columns=COLUMNS)
    locus_of = {name: locus for name, locus, _ in primers}
    for gid, index in indexes.items():
        for contig in index.seqs:
            sites = [s for name, _, seq in primers for s in find_sites(contig, index, seq, name, binding)]
            plus = sorted((s for s in sites if s.strand == "+"), key=lambda s: s.start)
            minus = sorted((s for s in sites if s.strand == "-"), key=lambda s: s.end)
            minus_ends = [s.end for s in minus]
            for p in plus:
                lo = bisect_right(minus_ends, p.end)
                hi = bisect_right(minus_ends, p.start + max_product_size)
                for m in minus[lo:hi]:
                    if m.start < p.end or locus_of[p.primer] == locus_of[m.primer]:
                        continue
                    rows.append({
                        "Genome": gid, "Contig": contig, "Start": p.start, "End": m.end, "Size": m.end - p.start,
                        "Forward_Primer": p.primer, "Reverse_Primer": m.primer,
                        "Forward_Locus": locus_of[p.primer], "Reverse_Locus": locus_of[m.primer],
                        "Forward_Tm": None if p.tm is None else round(p.tm, 1),
                        "Reverse_Tm": None if m.tm is None else round(m.tm, 1),
                    })
    return pd.DataFrame(rows, columns=COLUMNS)


def conflicting_loci(products: pd.DataFrame) -> List[Tuple[str, str]]:
    """Sorted, de-duplicated locus pairs that form at least one cross-pair product."""
    if products.empty:
        return []
    return sorted({tuple(sorted((a, b))) for a, b in zip(products["Forward_Locus"], products["Reverse_Locus"])})
