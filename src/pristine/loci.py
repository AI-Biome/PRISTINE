"""Locus discovery without annotation.

A target genome is tiled into overlapping windows of amplicon size. Each window is
mapped against every genome with minimap2. A window is kept as a candidate locus if,
in every genome, it has exactly one good hit and no substantial secondary hit
(single copy, present everywhere). The pipeline builds the locus alignments from the
``cs`` tags of these hits (see ``pseudoalign``); ``extract_loci`` writes the hit
regions as unaligned FASTA for use outside the pipeline.
"""

from __future__ import annotations

import logging
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .config import LocusDiscoveryParams
from .manifest import Genome
from .seqio import read_fasta_dict, revcomp, write_fasta

log = logging.getLogger(__name__)


@dataclass(frozen=True)
class Window:
    window_id: str
    contig: str
    start: int  # 0-based, inclusive
    end: int  # 0-based, exclusive


@dataclass
class PafHit:
    qname: str
    qlen: int
    qstart: int
    qend: int
    strand: str
    tname: str
    tlen: int
    tstart: int
    tend: int
    nmatch: int
    alnlen: int
    mapq: int
    primary: bool = True
    divergence: Optional[float] = None
    cs: str = ""

    @property
    def query_coverage(self) -> float:
        return (self.qend - self.qstart) / self.qlen if self.qlen else 0.0

    @property
    def identity(self) -> float:
        return self.nmatch / self.alnlen if self.alnlen else 0.0


@dataclass
class LocusHit:
    """Where a window lands in one genome."""

    genome_id: str
    contig: str
    start: int
    end: int
    strand: str
    divergence: float
    qstart: int = 0
    qend: int = 0
    cs: str = ""


@dataclass
class CandidateLocus:
    window: Window
    hits: Dict[str, LocusHit] = field(default_factory=dict)  # genome_id -> hit
    mean_nontarget_divergence: float = 0.0


def tile_genome(genome: Dict[str, str], window_size: int, step: int, max_n_fraction: float) -> List[Window]:
    windows: List[Window] = []
    for contig, seq in genome.items():
        n = len(seq)
        if n < window_size:
            continue
        starts = list(range(0, n - window_size + 1, step))
        if starts and starts[-1] + window_size < n:
            starts.append(n - window_size)  # cover the contig tail
        for s in starts:
            sub = seq[s : s + window_size]
            if sub.count("N") / window_size > max_n_fraction:
                continue
            windows.append(Window(f"{contig}:{s}-{s + window_size}", contig, s, s + window_size))
    return windows


def write_windows(windows: Sequence[Window], genome: Dict[str, str], path: Path) -> None:
    write_fasta(path, ((w.window_id, genome[w.contig][w.start : w.end]) for w in windows))


def minimap2_available() -> bool:
    return shutil.which("minimap2") is not None


def run_minimap2(windows_fa: Path, genome_fa: Path, paf_out: Path, threads: int, preset: str) -> None:
    if not minimap2_available():
        raise RuntimeError("minimap2 not found in PATH (run inside the pixi environment)")
    cmd = [
        "minimap2", "-c", "--cs", "-x", preset, "--secondary=yes", "-N", "10", "-p", "0.5",
        "-t", str(max(1, threads)), str(genome_fa), str(windows_fa),
    ]
    with open(paf_out, "w") as out, open(str(paf_out) + ".log", "w") as err:
        subprocess.run(cmd, stdout=out, stderr=err, check=True)


def parse_paf(path: Path) -> Dict[str, List[PafHit]]:
    hits: Dict[str, List[PafHit]] = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            hit = PafHit(
                qname=f[0], qlen=int(f[1]), qstart=int(f[2]), qend=int(f[3]), strand=f[4],
                tname=f[5], tlen=int(f[6]), tstart=int(f[7]), tend=int(f[8]),
                nmatch=int(f[9]), alnlen=int(f[10]), mapq=int(f[11]),
            )
            for tag in f[12:]:
                if tag.startswith("tp:A:"):
                    hit.primary = tag[5] == "P"
                elif tag.startswith("de:f:"):
                    hit.divergence = float(tag[5:])
                elif tag.startswith("cs:Z:"):
                    hit.cs = tag[5:]
            if hit.divergence is None:
                hit.divergence = 1.0 - hit.identity
            hits.setdefault(hit.qname, []).append(hit)
    return hits


def classify_hits(hits: List[PafHit], params: LocusDiscoveryParams) -> Optional[LocusHit]:
    """Return the single good hit, or None if the window is absent or multi-copy."""
    good = [h for h in hits if h.query_coverage >= params.min_query_coverage and h.divergence <= params.max_divergence]
    if len(good) != 1:
        return None
    best = good[0]
    others = [h for h in hits if h is not best and h.query_coverage >= params.max_secondary_coverage]
    if others:
        return None
    return LocusHit(genome_id="", contig=best.tname, start=best.tstart, end=best.tend,
                    strand=best.strand, divergence=best.divergence or 0.0,
                    qstart=best.qstart, qend=best.qend, cs=best.cs)


def select_candidate_loci(
    windows: Sequence[Window],
    paf_by_genome: Dict[str, Dict[str, List[PafHit]]],
    genomes: Sequence[Genome],
    params: LocusDiscoveryParams,
) -> List[CandidateLocus]:
    """Keep windows that are single-copy in every genome (non-targets optional)."""
    targets = [g for g in genomes if g.is_target]
    nontargets = [g for g in genomes if not g.is_target]
    kept: List[CandidateLocus] = []
    for w in windows:
        cand = CandidateLocus(window=w)
        ok = True
        for g in genomes:
            hit = classify_hits(paf_by_genome[g.genome_id].get(w.window_id, []), params)
            if hit is None:
                if g.is_target or params.require_in_nontargets:
                    ok = False
                    break
                continue
            hit.genome_id = g.genome_id
            cand.hits[g.genome_id] = hit
        if not ok or len(cand.hits) < len(targets):
            continue
        nt_div = [cand.hits[g.genome_id].divergence for g in nontargets if g.genome_id in cand.hits]
        cand.mean_nontarget_divergence = sum(nt_div) / len(nt_div) if nt_div else 0.0
        kept.append(cand)
    # Prefer loci that differ from the non-targets; cap the number to align.
    kept.sort(key=lambda c: -c.mean_nontarget_divergence)
    if len(kept) > params.max_candidate_loci:
        log.info("Capping candidate loci from %d to %d", len(kept), params.max_candidate_loci)
        kept = kept[: params.max_candidate_loci]
    return kept


def sanitize_contig(name: str, max_len: int = 40) -> str:
    """File-system safe, bounded-length version of a contig name."""
    clean = re.sub(r"[^A-Za-z0-9.]+", "_", name).strip("_")
    return clean[:max_len] or "contig"


def locus_id(window: Window) -> str:
    """Locus name used for files and tables: <contig>_<start>_<end>."""
    return f"{sanitize_contig(window.contig)}_{window.start}_{window.end}"


def locus_filename(window_id: str) -> str:
    return window_id.replace(":", "_").replace("-", "_") + ".fasta"


def extract_loci(
    candidates: Sequence[CandidateLocus],
    genome_seqs: Dict[str, Dict[str, str]],
    out_dir: Path,
) -> List[Tuple[CandidateLocus, Path]]:
    """Write one unaligned FASTA per locus with one record per genome (record id = genome_id)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for cand in candidates:
        records = []
        for gid, hit in cand.hits.items():
            seq = genome_seqs[gid][hit.contig][hit.start : hit.end]
            if hit.strand == "-":
                seq = revcomp(seq)
            records.append((gid, seq))
        path = out_dir / locus_filename(cand.window.window_id)
        write_fasta(path, records)
        written.append((cand, path))
    return written


def discover_loci(
    genomes: Sequence[Genome],
    reference: Genome,
    work_dir: Path,
    params: LocusDiscoveryParams,
    threads: int,
) -> Tuple[List[CandidateLocus], Dict[str, Dict[str, str]]]:
    """Full discovery for one group. Returns candidates and the loaded genome sequences."""
    work_dir.mkdir(parents=True, exist_ok=True)
    genome_seqs = {g.genome_id: read_fasta_dict(g.path) for g in genomes}
    windows = tile_genome(genome_seqs[reference.genome_id], params.window_size, params.step, params.max_n_fraction)
    log.info("Tiled reference %s into %d windows", reference.genome_id, len(windows))
    windows_fa = work_dir / "reference_windows.fasta"
    write_windows(windows, genome_seqs[reference.genome_id], windows_fa)

    # Several small minimap2 jobs at once use the cores better than one job with many
    # threads: indexing a bacterial genome is single-threaded, and so is the PAF parsing.
    jobs = max(1, min(len(genomes), threads // 2))
    per_job = max(1, threads // jobs)

    def map_one(g: Genome) -> Tuple[str, Dict[str, List[PafHit]]]:
        paf = work_dir / f"{g.genome_id}.paf"
        run_minimap2(windows_fa, g.path, paf, per_job, params.minimap2_preset)
        return g.genome_id, parse_paf(paf)

    with ThreadPoolExecutor(max_workers=jobs) as pool:
        paf_by_genome: Dict[str, Dict[str, List[PafHit]]] = dict(pool.map(map_one, genomes))
    candidates = select_candidate_loci(windows, paf_by_genome, genomes, params)
    log.info("%d of %d windows are single-copy and present in the required genomes", len(candidates), len(windows))
    return candidates, genome_seqs
