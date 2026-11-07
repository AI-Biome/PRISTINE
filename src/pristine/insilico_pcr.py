"""In silico PCR with a thermodynamic binding model.

Binding sites of a primer are found in two stages:

1. **Seed scan** (numpy): every ``seed_length``-mer inside the 3'-terminal
   ``three_prime_region`` bases of the primer is looked up in a cached k-mer index, so
   any site with at most one mismatch in that region is found (pigeonhole). Candidate
   sites are then kept if they have at most ``three_prime_max_mismatches`` mismatches in
   the 3' region and at most ``max_mismatches`` over the whole primer; the 5' part is
   also compared at a shift of +-1 base so that a template bulge in the 5' part does not
   hide a site (a bulge inside the 3' region is not modelled: such a 3' end is not
   extended efficiently). This is a cheap, deliberately permissive prefilter.
2. **Duplex Tm** (primer3 ``calc_heterodimer``): for every prefiltered site the melting
   temperature of the primer hybridised to the genomic window (with 3 extra template
   bases beyond the primer 5' end, so bulges and dangling ends are modelled) is
   computed. The site binds if that Tm is at least ``min_tm`` (annealing temperature
   minus a safety margin). With ``min_tm=None`` the mismatch count alone decides.

Products are all pairs of a plus-strand site and a downstream minus-strand site (any
primer combination, including F/F and R/R) up to ``max_product_size``. Every product
that does not lie inside the intended locus is an off-target product.
"""

from __future__ import annotations

from bisect import bisect_right
from dataclasses import dataclass, field
from typing import Dict, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np

from .seqio import revcomp

TEMPLATE_SLACK = 3  # extra template bases beyond the primer 5' end given to the thermodynamic model


class Thermo:
    """Thin wrapper around primer3's ThermoAnalysis with PCR-like salt conditions."""

    def __init__(self, mv_conc: float = 50.0, dv_conc: float = 1.5, dntp_conc: float = 0.6, dna_conc: float = 50.0):
        from primer3 import thermoanalysis

        self._ta = thermoanalysis.ThermoAnalysis(mv_conc=mv_conc, dv_conc=dv_conc, dntp_conc=dntp_conc, dna_conc=dna_conc)
        self._cache: Dict[Tuple[str, str], float] = {}

    def duplex_tm(self, primer: str, template_strand: str) -> float:
        """Tm of `primer` hybridised to `template_strand` (both written 5'->3'). Cached: the same
        genomic window recurs across related genomes."""
        key = (primer, template_strand)
        tm = self._cache.get(key)
        if tm is None:
            tm = float(self._ta.calc_heterodimer(primer, template_strand).tm)
            if len(self._cache) > 500_000:
                self._cache.clear()
            self._cache[key] = tm
        return tm

    def primer_tm(self, primer: str) -> float:
        return float(self._ta.calc_tm(primer))


@dataclass
class BindingParams:
    max_mismatches: int = 6  # over the whole primer
    three_prime_region: int = 8  # 3'-terminal bases scrutinised separately
    three_prime_max_mismatches: int = 1  # tolerated inside that region
    seed_length: int = 4  # exact k-mer used to find candidate sites
    min_tm: Optional[float] = None  # None: mismatch counts alone decide
    thermo: Optional[Thermo] = None

    def __post_init__(self) -> None:
        if self.min_tm is not None and self.thermo is None:
            self.thermo = Thermo()
        if self.seed_length < 1 or self.seed_length > 12:
            raise ValueError("seed_length must be between 1 and 12")


_BASE_CODE = np.full(256, 4, dtype=np.uint8)  # A,C,G,T -> 0..3; anything else -> 4 (invalid)
for _i, _b in enumerate(b"ACGT"):
    _BASE_CODE[_b] = _i


def kmer_codes(arr: np.ndarray, k: int) -> np.ndarray:
    """Integer code of every k-mer (2 bits per base); -1 where the k-mer contains a non-ACGT base."""
    n = arr.size
    if n < k:
        return np.empty(0, dtype=np.int32)
    codes = _BASE_CODE[arr].astype(np.int32)
    valid = codes < 4
    m = n - k + 1
    val = np.zeros(m, dtype=np.int32)
    ok = np.ones(m, dtype=bool)
    for j in range(k):
        val = val * 4 + codes[j : j + m]
        ok &= valid[j : j + m]
    val[~ok] = -1
    return val


def kmer_code(kmer: str) -> int:
    code = 0
    for ch in kmer.upper():
        v = int(_BASE_CODE[ord(ch)])
        if v == 4:
            return -2  # ambiguous primer base: never matches
        code = code * 4 + v
    return code


class GenomeIndex:
    """Upper-case sequences, their uint8 views and cached k-mer codes, built once per genome."""

    def __init__(self, seqs: Mapping[str, str]):
        self.seqs: Dict[str, str] = {c: s.upper() for c, s in seqs.items()}
        self.arrays: Dict[str, np.ndarray] = {c: np.frombuffer(s.encode("ascii"), dtype=np.uint8) for c, s in self.seqs.items()}
        self._kmers: Dict[Tuple[str, int], np.ndarray] = {}
        self._lookup: Dict[Tuple[str, int], Tuple[np.ndarray, np.ndarray]] = {}

    def kmers(self, contig: str, k: int) -> np.ndarray:
        key = (contig, k)
        if key not in self._kmers:
            self._kmers[key] = kmer_codes(self.arrays[contig], k)
        return self._kmers[key]

    def positions_of(self, contig: str, k: int, code: int) -> np.ndarray:
        key = (contig, k)
        if key not in self._lookup:
            codes = self.kmers(contig, k)
            order = np.argsort(codes, kind="stable").astype(np.int32)
            self._lookup[key] = (codes[order], order)
            self._kmers.pop(key, None)  # the sorted copy is all we need; keep memory per genome low
        sorted_codes, order = self._lookup[key]
        code32 = np.int32(code)  # same dtype as the array: a Python int would force a full int64 copy
        lo = np.searchsorted(sorted_codes, code32, side="left")
        hi = np.searchsorted(sorted_codes, code32, side="right")
        return order[lo:hi].astype(np.int64)

    def warm(self, k: int) -> "GenomeIndex":
        """Build the k-mer lookup of every contig now (e.g. before forking worker processes,
        so the workers share it copy-on-write instead of each building its own)."""
        for contig in self.seqs:
            self.positions_of(contig, k, 0)
        return self

    def __getitem__(self, contig: str) -> str:
        return self.seqs[contig]

    def items(self):
        return self.seqs.items()


@dataclass(frozen=True)
class Site:
    contig: str
    start: int  # 0-based inclusive
    end: int  # exclusive
    strand: str  # "+" primer reads 5'->3' along the contig, "-" along the reverse strand
    primer: str  # name ("F" / "R")
    mismatches: int
    tm: Optional[float] = None


@dataclass(frozen=True)
class Product:
    contig: str
    start: int
    end: int  # exclusive
    fwd_primer: str
    rev_primer: str
    fwd_mismatches: int
    rev_mismatches: int
    fwd_tm: Optional[float] = None
    rev_tm: Optional[float] = None

    @property
    def size(self) -> int:
        return self.end - self.start


def _as_bytes(s: str) -> np.ndarray:
    return np.frombuffer(s.encode("ascii"), dtype=np.uint8)


def _candidate_starts(index: GenomeIndex, contig: str, pattern: str, seed_offsets: Sequence[int], k: int) -> np.ndarray:
    """Start positions of windows where `pattern` has an exact k-mer at one of the seed offsets."""
    n = index.arrays[contig].size
    L = len(pattern)
    found = []
    for j in seed_offsets:
        code = kmer_code(pattern[j : j + k])
        if code < 0:
            continue
        found.append(index.positions_of(contig, k, code) - j)
    if not found:
        return np.empty(0, dtype=np.int64)
    starts = np.unique(np.concatenate(found))
    return starts[(starts >= 0) & (starts + L <= n)]


def _prefilter(
    arr: np.ndarray, starts: np.ndarray, pattern: np.ndarray, three_prime_cols: slice, head_cols: slice,
    max_mm: int, max_mm_3p: int, shifts: Sequence[int],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Cheap mismatch prefilter.

    The 3' region is compared in place (shift 0) and may hold at most `max_mm_3p`
    mismatches. The 5' part is compared at every shift in `shifts` (a template bulge in
    the 5' part shows up as a good match at a shift of one base) and the best count is
    used. Returns (kept starts, total mismatches, 3' mismatches)."""
    n = arr.size
    L = pattern.size
    empty = np.empty(0, dtype=np.int64)
    if starts.size == 0:
        return starts, empty, empty
    tp_off = np.arange(three_prime_cols.start, three_prime_cols.stop)
    tp = (arr[starts[:, None] + tp_off[None, :]] != pattern[three_prime_cols][None, :]).sum(axis=1)
    keep = tp <= max_mm_3p
    starts, tp = starts[keep], tp[keep]
    if starts.size == 0:
        return starts, empty, empty
    head_off = np.arange(head_cols.start, head_cols.stop)
    head_best = np.zeros(starts.size, dtype=np.int64)
    if head_off.size:
        head_best = np.full(starts.size, L + 1, dtype=np.int64)
        for d in shifts:
            s = starts + d
            ok = (s + head_cols.start >= 0) & (s + head_cols.stop <= n)
            if not ok.any():
                continue
            mm = np.full(starts.size, L + 1, dtype=np.int64)
            mm[ok] = (arr[s[ok][:, None] + head_off[None, :]] != pattern[head_cols][None, :]).sum(axis=1)
            head_best = np.minimum(head_best, mm)
    total = head_best + tp
    keep = total <= max_mm
    return starts[keep], total[keep], tp[keep]


def find_sites(contig: str, genome: Union[str, GenomeIndex], primer: str, name: str, params: BindingParams) -> List[Site]:
    """All binding sites of `primer` on both strands of one contig."""
    index = genome if isinstance(genome, GenomeIndex) else GenomeIndex({contig: genome})
    seq, arr = index.seqs[contig], index.arrays[contig]
    primer = primer.upper()
    L = len(primer)
    R = max(1, min(params.three_prime_region, L))
    k = max(1, min(params.seed_length, R))
    shifts = (0, -1, 1) if params.min_tm is not None else (0,)
    n = arr.size
    sites: List[Site] = []

    def accept(template_window: str) -> Optional[float]:
        if params.min_tm is None:
            return None
        tm = params.thermo.duplex_tm(primer, template_window)
        return tm if tm >= params.min_tm else -1.0

    # plus strand: primer reads along the contig; its 3' end is the last base of the window
    starts = _candidate_starts(index, contig, primer, range(L - R, L - k + 1), k)
    starts, tot, _tp = _prefilter(arr, starts, _as_bytes(primer), slice(L - R, L), slice(0, L - R),
                                  params.max_mismatches, params.three_prime_max_mismatches, shifts)
    for start, mm in zip(starts.tolist(), tot.tolist()):
        window = seq[max(0, start - TEMPLATE_SLACK) : start + L]
        tm = accept(revcomp(window))
        if tm is None or tm >= 0:
            sites.append(Site(contig, start, start + L, "+", name, mm, tm))

    # minus strand: the reverse complement reads along the contig; the primer 3' end is its first base
    rc = revcomp(primer)
    starts = _candidate_starts(index, contig, rc, range(0, R - k + 1), k)
    starts, tot, _tp = _prefilter(arr, starts, _as_bytes(rc), slice(0, R), slice(R, L),
                                  params.max_mismatches, params.three_prime_max_mismatches, shifts)
    for start, mm in zip(starts.tolist(), tot.tolist()):
        window = seq[start : min(n, start + L + TEMPLATE_SLACK)]
        tm = accept(window)
        if tm is None or tm >= 0:
            sites.append(Site(contig, start, start + L, "-", name, mm, tm))
    return sites


def amplify(
    genome: Union[Mapping[str, str], GenomeIndex],
    forward: str,
    reverse: str,
    max_product_size: int,
    params: Optional[BindingParams] = None,
) -> List[Product]:
    params = params or BindingParams()
    index = genome if isinstance(genome, GenomeIndex) else GenomeIndex(genome)
    products: List[Product] = []
    for contig in index.seqs:
        sites = find_sites(contig, index, forward, "F", params) + find_sites(contig, index, reverse, "R", params)
        plus = sorted((s for s in sites if s.strand == "+"), key=lambda s: s.start)
        minus = sorted((s for s in sites if s.strand == "-"), key=lambda s: s.end)
        minus_ends = [s.end for s in minus]
        for p in plus:
            lo = bisect_right(minus_ends, p.end)
            hi = bisect_right(minus_ends, p.start + max_product_size)
            for m in minus[lo:hi]:
                if m.start >= p.end:  # non-overlapping primers
                    products.append(Product(contig, p.start, m.end, p.primer, m.primer, p.mismatches, m.mismatches, p.tm, m.tm))
    return products


@dataclass
class GenomePcrResult:
    genome_id: str
    products: List[Product] = field(default_factory=list)
    on_target: List[Product] = field(default_factory=list)
    off_target: List[Product] = field(default_factory=list)

    @property
    def n_products(self) -> int:
        return len(self.products)


def classify_products(
    genome_id: str,
    products: Sequence[Product],
    locus_region: Optional[Tuple[str, int, int]],
    margin: int = 50,
) -> GenomePcrResult:
    """Split products into on-target (contained in the locus region) and off-target."""
    on, off = [], []
    for pr in products:
        if locus_region is not None:
            contig, lo, hi = locus_region
            if pr.contig == contig and pr.start >= lo - margin and pr.end <= hi + margin:
                on.append(pr)
                continue
        off.append(pr)
    return GenomePcrResult(genome_id, list(products), on, off)


def product_sequence(genome: Union[Mapping[str, str], GenomeIndex], product: Product) -> str:
    return genome[product.contig][product.start : product.end]
