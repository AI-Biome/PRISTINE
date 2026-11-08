"""Informative-site scoring of a locus alignment.

A column is *informative* for non-target species S if the set of states shown by the
target group (all target genomes) is disjoint from the set of states shown by the
genomes of S. Ambiguous bases are ignored. States are bit flags, so a whole alignment
is scored with a few numpy operations. A gap is a state, but a run of
consecutive informative gap columns is one *indel event* and is counted once (at its
first column); otherwise a single 200 bp insertion would look like 200 SNPs. A 1-nt
indel inside a homopolymer of >= HOMOPOLYMER_WEAK bases (in the allele with the base)
is a *weak* indel: the typical nanopore error, so reads cannot be told apart by it. Weak
indels are reported but not counted as informative sites. Internal
informativeness (between target species of the same group) is scored the same way
pairwise, so groups with several target species also get species-level markers.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Mapping, Sequence, Set, Tuple

import numpy as np

from .manifest import Genome

VALID_STATES = frozenset("ACGT-")
HOMOPOLYMER_WEAK = 4  # 1-nt indels in homopolymers this long are not counted


@dataclass
class SiteSet:
    """Informative sites between two groups of sequences."""

    substitutions: List[int] = field(default_factory=list)  # alignment columns
    indel_events: List[int] = field(default_factory=list)  # first column of each gap run
    weak_indels: List[int] = field(default_factory=list)  # 1-nt indels in homopolymers, not counted

    @property
    def positions(self) -> List[int]:
        return sorted(set(self.substitutions) | set(self.indel_events))

    @property
    def count(self) -> int:
        return len(self.substitutions) + len(self.indel_events)


@dataclass
class LocusScore:
    locus: str
    aln_len: int
    n_target_genomes: int
    n_nontarget_genomes: int
    external: Dict[str, SiteSet] = field(default_factory=dict)  # non-target species -> sites
    internal: Dict[Tuple[str, str], SiteSet] = field(default_factory=dict)  # (target sp, other sp)
    med_target_len: float = 0.0
    med_nontarget_len: float = 0.0

    # -- convenience views ------------------------------------------------------
    @property
    def external_positions(self) -> Dict[str, List[int]]:
        return {sp: s.positions for sp, s in self.external.items()}

    @property
    def internal_positions(self) -> Dict[Tuple[str, str], List[int]]:
        return {k: s.positions for k, s in self.internal.items()}

    def _prop(self, s: SiteSet) -> float:
        return s.count / self.aln_len if self.aln_len else 0.0

    @property
    def avg_prop_external(self) -> float:
        props = [self._prop(s) for s in self.external.values()]
        return float(np.mean(props)) if props else 0.0

    @property
    def avg_prop_internal(self) -> float:
        if not self.internal:
            return 0.0
        per_target: Dict[str, List[float]] = {}
        for (t, _o), s in self.internal.items():
            per_target.setdefault(t, []).append(self._prop(s))
        return float(np.mean([np.mean(v) for v in per_target.values()]))

    @property
    def combined_score(self) -> float:
        if self.external and self.internal:
            return 0.5 * (self.avg_prop_external + self.avg_prop_internal)
        if self.external:
            return self.avg_prop_external
        return self.avg_prop_internal

    @property
    def min_external_count(self) -> int:
        return min((s.count for s in self.external.values()), default=0)

    def all_positions(self) -> List[int]:
        out: Set[int] = set()
        for s in list(self.external.values()) + list(self.internal.values()):
            out.update(s.positions)
        return sorted(out)

    def positions_by_key(self) -> Dict[str, List[int]]:
        out = {f"External_{sp}": s.positions for sp, s in self.external.items()}
        out.update({f"Internal_{t}_vs_{o}": s.positions for (t, o), s in self.internal.items()})
        return out

    def count_in_region(self, start: int, end: int) -> Dict[str, int]:
        """Informative sites inside [start, end] (alignment coords), per key."""
        return {k: sum(start <= x <= end for x in p) for k, p in self.positions_by_key().items()}

    def to_record(self) -> Dict[str, object]:
        rec: Dict[str, object] = {
            "Locus": self.locus,
            "Alignment_Length": self.aln_len,
            "Combined_Informativeness_Score": self.combined_score,
            "Avg_Prop_Informative_SNPs_External": self.avg_prop_external,
            "Avg_Prop_Informative_SNPs_Internal": self.avg_prop_internal,
            "Min_Abs_External": self.min_external_count,
            "Med_Target_Len": self.med_target_len,
            "Med_NonTarget_Len": self.med_nontarget_len,
            "Len_Diff": self.med_target_len - self.med_nontarget_len,
        }
        for sp, s in sorted(self.external.items()):
            self._add_site_columns(rec, f"External_{sp}", s)
        for (t, o), s in sorted(self.internal.items()):
            self._add_site_columns(rec, f"Internal_{t}_vs_{o}", s)
        return rec

    def _add_site_columns(self, rec: Dict[str, object], key: str, s: SiteSet) -> None:
        rec[f"Prop_{key}"] = self._prop(s)
        rec[f"Abs_{key}"] = s.count
        rec[f"Subst_{key}"] = len(s.substitutions)
        rec[f"Indel_{key}"] = len(s.indel_events)
        rec[f"Weak_Indel_{key}"] = len(s.weak_indels)
        rec[f"SNP_Pos_{key}"] = ",".join(map(str, s.positions))
        rec[f"Core_Region_{key}"] = core_region(s.positions, self.aln_len)


def core_region(positions: Sequence[int], aln_len: int) -> str:
    if len(positions) < 2:
        return ""
    mean, std = float(np.mean(positions)), float(np.std(positions))
    lower = max(0, int(mean - 2 * std))
    upper = min(aln_len - 1, int(mean + 2 * std))
    return f"{lower}-{upper}"


# Base states as bit flags; ambiguous bases are 0 and therefore ignored.
_STATE_BITS = np.zeros(256, dtype=np.uint8)
for _ch, _bit in ((b"A", 1), (b"C", 2), (b"G", 4), (b"T", 8), (b"-", 16), (b"a", 1), (b"c", 2), (b"g", 4), (b"t", 8)):
    _STATE_BITS[_ch[0]] = _bit
_GAP_BIT = 16


def _states(bases: Iterable[str]) -> Set[str]:
    return {b for b in bases if b in VALID_STATES}


def is_informative(group_states: Set[str], other_states: Set[str]) -> bool:
    return bool(group_states) and bool(other_states) and group_states.isdisjoint(other_states)


def _state_matrix(seqs: Sequence[str]) -> np.ndarray:
    return np.vstack([_STATE_BITS[np.frombuffer(s.encode("ascii", "replace"), dtype=np.uint8)] for s in seqs])


def _or_rows(mat: np.ndarray, rows: Sequence[int]) -> np.ndarray:
    return np.bitwise_or.reduce(mat[list(rows)], axis=0)


def homopolymer_run(seq: str, col: int) -> int:
    """Length of the homopolymer containing the base at alignment column `col` of `seq`,
    gaps skipped (0 if `seq` has a gap there)."""
    b = seq[col].upper()
    if b not in "ACGT":
        return 0
    n = 1
    for step in (-1, 1):
        i = col + step
        while 0 <= i < len(seq):
            c = seq[i].upper()
            if c == "-":
                i += step
                continue
            if c != b:
                break
            n += 1
            i += step
    return n


def sites_between(group_mask: np.ndarray, other_mask: np.ndarray, seqs: Sequence[str] = ()) -> SiteSet:
    """Informative columns between two groups given their per-column state bitmasks.

    Consecutive informative columns that involve a gap are one indel event. With the
    aligned `seqs`, a 1-column indel event whose base lies in a homopolymer of
    >= HOMOPOLYMER_WEAK in any sequence is weak and not counted."""
    informative = (group_mask != 0) & (other_mask != 0) & ((group_mask & other_mask) == 0)
    gap = ((group_mask | other_mask) & _GAP_BIT) != 0
    subst = informative & ~gap
    indel = informative & gap
    prev = np.concatenate(([False], indel[:-1]))
    nxt = np.concatenate((indel[1:], [False]))
    first_of_run = indel & ~prev
    events = np.flatnonzero(first_of_run).tolist()
    weak = []
    if len(seqs):
        single = set(np.flatnonzero(first_of_run & ~nxt).tolist())
        weak = [c for c in events if c in single and max(homopolymer_run(s, c) for s in seqs) >= HOMOPOLYMER_WEAK]
        events = [c for c in events if c not in set(weak)]
    return SiteSet(substitutions=np.flatnonzero(subst).tolist(), indel_events=events, weak_indels=weak)


def score_alignment(locus: str, aln: Mapping[str, str], genomes: Mapping[str, Genome]) -> LocusScore:
    """Score one alignment {genome_id: aligned_sequence} (vectorised over columns)."""
    ids = [gid for gid in aln if gid in genomes]
    if not ids:
        raise ValueError(f"{locus}: no alignment record matches a known genome id")
    seqs = [aln[gid] for gid in ids]
    aln_len = len(seqs[0])
    if any(len(s) != aln_len for s in seqs):
        raise ValueError(f"{locus}: sequences have unequal aligned lengths")
    row = {gid: i for i, gid in enumerate(ids)}
    mat = _state_matrix(seqs)

    target_ids = [gid for gid in ids if genomes[gid].is_target]
    nontarget_species: Dict[str, List[str]] = {}
    for gid in ids:
        if not genomes[gid].is_target:
            nontarget_species.setdefault(genomes[gid].species, []).append(gid)
    target_species: Dict[str, List[str]] = {}
    for gid in target_ids:
        target_species.setdefault(genomes[gid].species, []).append(gid)

    external: Dict[str, SiteSet] = {}
    if target_ids:
        group_mask = _or_rows(mat, [row[g] for g in target_ids])
        for sp, gids in nontarget_species.items():
            external[sp] = sites_between(group_mask, _or_rows(mat, [row[g] for g in gids]), seqs)
    internal: Dict[Tuple[str, str], SiteSet] = {}
    if len(target_species) > 1:
        sp_mask = {sp: _or_rows(mat, [row[g] for g in gids]) for sp, gids in target_species.items()}
        for t in target_species:
            for o in target_species:
                if t != o:
                    internal[(t, o)] = sites_between(sp_mask[t], sp_mask[o], seqs)

    ungapped = {gid: len(s) - s.count("-") for gid, s in zip(ids, seqs)}
    tl = [ungapped[g] for g in target_ids]
    nl = [ungapped[g] for g in ids if g not in target_ids]
    return LocusScore(
        locus=locus,
        aln_len=aln_len,
        n_target_genomes=len(target_ids),
        n_nontarget_genomes=len(ids) - len(target_ids),
        external=external,
        internal=internal,
        med_target_len=float(np.median(tl)) if tl else 0.0,
        med_nontarget_len=float(np.median(nl)) if nl else 0.0,
    )
