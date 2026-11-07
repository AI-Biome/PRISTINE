"""Primer design with primer3 on an N-masked template, forcing the amplicon to cover
the SNP-dense region of the locus."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import primer3

log = logging.getLogger(__name__)

GAP_TOLERANCE = 20  # bp; SNP windows closer than this are merged into one target


@dataclass
class PrimerPair:
    locus: str
    pair_index: int
    left_seq: str
    right_seq: str
    left_start: int  # template (alignment) coordinate of the left primer's 5' base
    left_len: int
    right_end: int  # template coordinate of the right primer's 5' base (= last base of amplicon)
    right_len: int
    left_tm: float
    right_tm: float
    left_gc: float
    right_gc: float
    template_product_size: int
    penalty: float
    snps_in_amplicon: Dict[str, int] = field(default_factory=dict)

    @property
    def amplicon_start(self) -> int:
        return self.left_start

    @property
    def amplicon_end(self) -> int:  # inclusive
        return self.right_end


def best_window(positions: Sequence[int], window_size: int) -> Tuple[Optional[int], int]:
    """Start of the window of `window_size` bp holding the most positions, and that count."""
    pos = sorted(set(positions))
    best_start, best_count = None, 0
    for i, start in enumerate(pos):
        end = start + window_size
        count = 0
        for x in pos[i:]:
            if x > end:
                break
            count += 1
        if count > best_count:
            best_start, best_count = start, count
    return best_start, best_count


def choose_target_region(
    positions_by_key: Mapping[str, Sequence[int]],
    window_size: int,
    min_snps: int,
    max_target_span: int,
    template_len: int,
) -> Optional[Tuple[int, int]]:
    """Return (start, length) of the region the amplicon must cover, or None.

    One SNP-dense window per non-target species (or target-species pair) is chosen;
    overlapping windows are merged. If the merged region is wider than
    `max_target_span`, the single densest window over all positions is used instead.
    """
    windows: List[Tuple[int, int]] = []
    for positions in positions_by_key.values():
        start, count = best_window(positions, window_size)
        if start is not None and count >= min_snps:
            windows.append((start, min(start + window_size, template_len)))
    if not windows:
        return None
    windows.sort()
    merged = [list(windows[0])]
    for s, e in windows[1:]:
        if s <= merged[-1][1] + GAP_TOLERANCE:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    span_start, span_end = merged[0][0], merged[-1][1]
    if span_end - span_start > max_target_span:
        all_pos = sorted({p for ps in positions_by_key.values() for p in ps})
        start, count = best_window(all_pos, max_target_span)
        if start is None or count < min_snps:
            return None
        span_start, span_end = start, min(start + max_target_span, template_len)
    return span_start, span_end - span_start


def design_primers(
    locus: str,
    template: str,
    target: Tuple[int, int],
    primer3_params: Mapping[str, Any],
) -> Tuple[List[PrimerPair], str]:
    """Run primer3 and return candidate pairs plus primer3's explanation string."""
    seq_args = {
        "SEQUENCE_ID": locus,
        "SEQUENCE_TEMPLATE": template,
        "SEQUENCE_TARGET": [list(target)],
    }
    result = primer3.bindings.design_primers(seq_args, dict(primer3_params))
    explain = "; ".join(
        f"{k}: {result.get(k)}" for k in ("PRIMER_LEFT_EXPLAIN", "PRIMER_RIGHT_EXPLAIN", "PRIMER_PAIR_EXPLAIN") if k in result
    )
    pairs: List[PrimerPair] = []
    for i in range(int(result.get("PRIMER_PAIR_NUM_RETURNED", 0))):
        left_start, left_len = result[f"PRIMER_LEFT_{i}"]
        right_end, right_len = result[f"PRIMER_RIGHT_{i}"]
        pairs.append(
            PrimerPair(
                locus=locus,
                pair_index=i,
                left_seq=result[f"PRIMER_LEFT_{i}_SEQUENCE"],
                right_seq=result[f"PRIMER_RIGHT_{i}_SEQUENCE"],
                left_start=int(left_start),
                left_len=int(left_len),
                right_end=int(right_end),
                right_len=int(right_len),
                left_tm=float(result[f"PRIMER_LEFT_{i}_TM"]),
                right_tm=float(result[f"PRIMER_RIGHT_{i}_TM"]),
                left_gc=float(result[f"PRIMER_LEFT_{i}_GC_PERCENT"]),
                right_gc=float(result[f"PRIMER_RIGHT_{i}_GC_PERCENT"]),
                template_product_size=int(result[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"]),
                penalty=float(result[f"PRIMER_PAIR_{i}_PENALTY"]),
            )
        )
    return pairs, explain
