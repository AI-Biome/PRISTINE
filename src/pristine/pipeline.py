"""Orchestration of the PRISTINE pipeline for every target group."""

from __future__ import annotations

import logging
import multiprocessing
import shutil
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import pandas as pd

from . import __version__
from .alignment import align_loci
from .config import Config
from .consensus import iupac_consensus, masked_template
from .insilico_pcr import (
    BindingParams,
    GenomeIndex,
    GenomePcrResult,
    Product,
    Thermo,
    amplify,
    classify_products,
    is_near_miss,
    product_sequence,
)
from .loci import CandidateLocus, Window, discover_loci, locus_id
from .manifest import Genome, choose_reference, discover_groups, is_fasta, species_from_stem
from .primers import PrimerPair, choose_target_region, design_primers
from .pseudoalign import build_alignment
from .scoring import HOMOPOLYMER_WEAK, LocusScore, score_alignment
from .seqio import read_fasta_dict, revcomp, write_fasta

log = logging.getLogger(__name__)


@dataclass
class PairEvaluation:
    pair: PrimerPair
    passed: bool
    reasons: List[str] = field(default_factory=list)
    per_genome: Dict[str, GenomePcrResult] = field(default_factory=dict)
    target_amplicons: List[str] = field(default_factory=list)  # set before the background screen
    background_products: List[str] = field(default_factory=list)  # accepted (distinguishable) ones
    background_sequences: Dict[str, str] = field(default_factory=dict)  # closest accepted product per genome


def list_background(background_dir: Optional[Path]) -> List[Path]:
    if background_dir is None:
        return []
    if not background_dir.is_dir():
        raise FileNotFoundError(f"primer_selection.background_dir not found: {background_dir}")
    return [p for p in sorted(background_dir.iterdir()) if is_fasta(p)]


def binding_params_from_config(cfg: Config) -> BindingParams:
    sel = cfg.primer_selection
    p3 = cfg.primer3_params
    thermo = None
    if cfg.min_site_tm is not None:
        thermo = Thermo(
            mv_conc=float(p3.get("PRIMER_SALT_MONOVALENT", 50.0)),
            dv_conc=float(p3.get("PRIMER_SALT_DIVALENT", 1.5)),
            dntp_conc=float(p3.get("PRIMER_DNTP_CONC", 0.6)),
            dna_conc=float(p3.get("PRIMER_DNA_CONC", 250.0)),
        )
    return BindingParams(max_mismatches=sel.max_mismatches, three_prime_region=sel.three_prime_region,
                         three_prime_max_mismatches=sel.three_prime_max_mismatches, seed_length=sel.seed_length,
                         min_tm=cfg.min_site_tm, thermo=thermo, near_miss_tm=cfg.near_miss_site_tm)


def split_near_misses(products: Sequence[Product], binding: BindingParams) -> Tuple[List[Product], List[Product]]:
    """(products that count, near misses) of an `amplify(..., include_near_misses=True)` result."""
    real, near = [], []
    for pr in products:
        (near if is_near_miss(pr, binding) else real).append(pr)
    return real, near


def describe_near_misses(products: Sequence[Product]) -> str:
    """Failure reason for near-miss products: count, sizes and the Tm of the weaker site."""
    sizes = ",".join(str(pr.size) for pr in sorted(products, key=lambda p: p.size))
    tms = [pr.weakest_tm for pr in products if pr.weakest_tm is not None]
    text = f"{len(products)} near-miss product(s) of {sizes} bp"
    if tms:
        text += f", weaker site Tm up to {max(tms):.1f} C"
    return text


# ---- parallel in silico PCR ---------------------------------------------------------
# Worker processes are forked after the parent has stored the genome indexes (with their
# k-mer lookups already built) in _SHARED, so every worker reads the same memory
# copy-on-write instead of receiving or rebuilding multi-megabyte indexes. The work is
# split by primer pair: the duplex-Tm cache is keyed by primer, so a worker that runs one
# pair over all genomes gets the full benefit of the cache.
_SHARED: Dict[str, object] = {}


def _fork_pool_possible(workers: int) -> bool:
    # fork is the only start method that shares the indexes without copying them; it is
    # the default and safe on Linux, but not on macOS, where PRISTINE runs serially.
    return workers >= 2 and sys.platform.startswith("linux")


def _fork_pool(workers: int) -> Optional[ProcessPoolExecutor]:
    """A fork-based process pool, or None (run serially) for one worker or off Linux."""
    if not _fork_pool_possible(workers):
        return None
    return ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context("fork"))


def _seed_k(binding: BindingParams) -> int:
    """k-mer size of the seed scan, as find_sites derives it for primers of normal length."""
    return max(1, min(binding.seed_length, binding.three_prime_region))


def _evaluate_job(job: Tuple[PrimerPair, str]) -> "PairEvaluation":
    pair, locus = job
    s = _SHARED
    return evaluate_pair(pair, s["genomes"], s["indexes"], s["cands"][locus], s["cfg"], s["binding"])


def _background_job(job: Tuple[str, str, Sequence[str], bool]) -> "_BackgroundResult":
    left, right, amplicons, distinguish = job
    s = _SHARED
    return _screen_one(s["bg_index"], left, right, amplicons, distinguish, s["max_product"], s["binding"], s["cap"])


@dataclass
class _BackgroundResult:
    products: List[Product]  # products that count
    diff: Optional[int] = None  # fewest differences of `products` from a target amplicon
    closest: Optional[str] = None  # sequence of that product
    near: List[Product] = field(default_factory=list)  # near misses: reject the pair


def _closest(index: GenomeIndex, products: Sequence[Product], amplicons: Sequence[str], cap: int) -> Tuple[int, str]:
    diffs = [(min_differences(seq, amplicons, cap), seq)
             for seq in (product_sequence(index, pr) for pr in products)]
    return min(diffs, key=lambda d: d[0])


def _screen_one(index: GenomeIndex, left: str, right: str, amplicons: Sequence[str], distinguish: bool,
                max_product_size: int, binding: BindingParams, cap: int) -> _BackgroundResult:
    """Products and near misses of one pair in one background genome and, with
    `distinguish`, the fewest differences between the products and the target amplicons
    (not needed when a near miss already rejects the pair)."""
    products, near = split_near_misses(amplify(index, left, right, max_product_size, binding,
                                               include_near_misses=True), binding)
    res = _BackgroundResult(products, near=near)
    if distinguish and products and not near:
        res.diff, res.closest = _closest(index, products, amplicons, cap)
    return res


# ---- distinguishing background products by sequence ---------------------------------
_KMER = 15
MANY = 10**6  # "differs at many sites"; returned without aligning
_ALIGNER = None


def _aligner():
    global _ALIGNER
    if _ALIGNER is None:
        from Bio.Align import PairwiseAligner
        _ALIGNER = PairwiseAligner(mode="global", match_score=1, mismatch_score=-2,
                                   open_gap_score=-4, extend_gap_score=-1)
    return _ALIGNER


def _hp_run(seq: str, start: int, end: int) -> int:
    """Homopolymer length around seq[start:end] if that is a single repeated base, else 1."""
    if end - start != 1:
        return 1
    b = seq[start]
    i, j = start, end
    while i > 0 and seq[i - 1] == b:
        i -= 1
    while j < len(seq) and seq[j] == b:
        j += 1
    return j - i


def _events(a: str, b: str) -> int:
    """Substitutions plus indel events (a gap run counts once) in a global alignment.
    A 1-nt indel in a homopolymer of >= HOMOPOLYMER_WEAK is not counted (nanopore error)."""
    aln = _aligner().align(a, b)[0]
    ta, qa = aln.aligned
    if not len(ta):
        return 1
    n = sum(x != y for (ts, te), (qs, qe) in zip(ta, qa) for x, y in zip(a[ts:te], b[qs:qe]))
    # gaps: before the first block, between blocks, after the last block
    bounds = [((0, ta[0][0]), (0, qa[0][0]))]
    bounds += [((ta[k - 1][1], ta[k][0]), (qa[k - 1][1], qa[k][0])) for k in range(1, len(ta))]
    bounds.append(((ta[-1][1], len(a)), (qa[-1][1], len(b))))
    for (as_, ae), (bs, be) in bounds:
        if ae == as_ and be == bs:
            continue
        weak = (be == bs and _hp_run(a, as_, ae) >= HOMOPOLYMER_WEAK) or \
               (ae == as_ and _hp_run(b, bs, be) >= HOMOPOLYMER_WEAK)
        n += not weak
    return n


def min_differences(product: str, amplicons: Sequence[str], cap: int) -> int:
    """Fewest differences (substitutions + indel events) between `product`, in either
    orientation, and any of the target `amplicons`. `cap` is the decision threshold
    (min_snps); the k-mer shortcut below only skips sequences far beyond it.

    Sequences sharing fewer than half of their 15-mers are not aligned and reported as
    MANY: fewer than `cap` events could only make them that different through indels of
    hundreds of bases, which reads tell apart anyway."""
    best = MANY
    for seq in (product.upper(), revcomp(product.upper())):
        kmers = {seq[i : i + _KMER] for i in range(len(seq) - _KMER + 1)}
        for amp in amplicons:
            if seq == amp:
                return 0
            amp_kmers = {amp[i : i + _KMER] for i in range(len(amp) - _KMER + 1)}
            if not kmers or len(kmers & amp_kmers) < 0.5 * min(len(kmers), len(amp_kmers)):
                continue
            best = min(best, _events(seq, amp))
    return best


def evaluate_pair(
    pair: PrimerPair,
    genomes: Sequence[Genome],
    indexes: Mapping[str, GenomeIndex],
    locus: CandidateLocus,
    cfg: Config,
    binding: BindingParams,
) -> PairEvaluation:
    """Informative-site criterion, then in silico PCR of one pair against the group's genomes.

    The site criterion needs no PCR, so a pair that fails it is not amplified at all (its
    product sizes stay empty in the reports)."""
    sel = cfg.primer_selection
    ev = PairEvaluation(pair=pair, passed=True)
    size_lo = int(cfg.min_product_size * 0.9)
    size_hi = int(cfg.max_product_size * 1.1)

    external = {k: v for k, v in pair.snps_in_amplicon.items() if k.startswith("External_")}
    min_snps = cfg.snp_primer_design.min_snps
    if external:
        lacking = [k for k, v in external.items() if v < min_snps]
        if cfg.snp_primer_design.require_snps_for_all_nontargets and lacking:
            ev.passed = False
            ev.reasons.append("amplicon has < %d informative sites for %s" % (min_snps, ", ".join(lacking)))
        elif not cfg.snp_primer_design.require_snps_for_all_nontargets and len(lacking) == len(external):
            ev.passed = False
            ev.reasons.append(f"amplicon has < {min_snps} informative sites for every non-target")
    elif sum(pair.snps_in_amplicon.values()) < min_snps:
        ev.passed = False
        ev.reasons.append(f"amplicon has < {min_snps} informative sites")
    if not ev.passed:
        ev.reasons.append("in silico PCR skipped")
        return ev

    for g in genomes:
        hit = locus.hits.get(g.genome_id)
        region = (hit.contig, hit.start, hit.end) if hit else None
        products, near = split_near_misses(amplify(indexes[g.genome_id], pair.left_seq, pair.right_seq,
                                                   sel.max_offtarget_product_size, binding, include_near_misses=True),
                                            binding)
        near_off = classify_products(g.genome_id, near, region).off_target
        if near_off:
            ev.passed = False
            ev.reasons.append(f"{g.genome_id}: {describe_near_misses(near_off)}")
        res = classify_products(g.genome_id, products, region)
        ev.per_genome[g.genome_id] = res
        if res.off_target:
            ev.passed = False
            ev.reasons.append(f"{g.genome_id}: {len(res.off_target)} off-target product(s)")
        n_on = len(res.on_target)
        # A non-target without the locus (require_in_nontargets: false) must give no product
        # at all; any product it gives is off-target and already failed above.
        if (g.is_target or (sel.conservation == "all" and hit is not None)) and n_on == 0:
            ev.passed = False
            ev.reasons.append(f"{g.genome_id}: no product at locus")
        if n_on > 1:
            ev.passed = False
            ev.reasons.append(f"{g.genome_id}: {n_on} products at locus")
        for pr in res.on_target:
            if not (size_lo <= pr.size <= size_hi):
                ev.passed = False
                ev.reasons.append(f"{g.genome_id}: product size {pr.size} outside {size_lo}-{size_hi}")
    return ev


def screen_background(
    evaluations: Sequence[PairEvaluation],
    background_files: Sequence[Path],
    binding: BindingParams,
    max_product_size: int,
    workers: int = 1,
    mode: str = "reject",
    min_snps: int = 3,
) -> None:
    """One streaming pass over the background genomes.

    With mode "reject" any background product fails the pair. With "distinguish" a
    product only fails it if it differs from every target amplicon (`ev.target_amplicons`)
    by fewer than `min_snps` substitutions + indel events, i.e. if reads of it could be
    mistaken for the target. Genomes are loaded one at a time; within a genome the
    still-passing pairs are split over `workers` processes that share the genome's index."""
    distinguish = mode == "distinguish"
    alive = [ev for ev in evaluations if ev.passed]
    log.info("Screening %d candidate pair(s) against %d background genome(s)", len(alive), len(background_files))
    t_start = time.time()
    for i, path in enumerate(background_files, 1):
        alive = [ev for ev in alive if ev.passed]
        if not alive:
            break
        index = GenomeIndex(read_fasta_dict(path))
        bid = path.stem
        jobs = [(ev.pair.left_seq, ev.pair.right_seq, ev.target_amplicons, distinguish) for ev in alive]
        pool = _fork_pool(min(workers, len(alive)))
        if pool is None:
            results = [_screen_one(index, l, r, a, d, max_product_size, binding, min_snps) for l, r, a, d in jobs]
        else:
            index.warm(_seed_k(binding))  # build once here, not once per worker
            _SHARED.update(bg_index=index, binding=binding, max_product=max_product_size, cap=min_snps)
            try:
                with pool:
                    results = list(pool.map(_background_job, jobs))
            finally:
                _SHARED.pop("bg_index", None)
        for ev, res in zip(alive, results):
            products, diff, closest = res.products, res.diff, res.closest
            if res.near:
                ev.passed = False
                ev.per_genome[bid] = classify_products(bid, res.near, None)
                ev.reasons.append(f"background {bid}: {describe_near_misses(res.near)}")
                continue
            if not products:
                continue
            if distinguish and diff >= min_snps:
                ev.background_products.append(
                    f"{bid}:{len(products)} ({'many' if diff >= MANY else diff} differences)")
                ev.background_sequences[bid] = closest
                continue
            ev.passed = False
            ev.per_genome[bid] = classify_products(bid, products, None)
            if distinguish:
                ev.reasons.append(f"background {bid}: {len(products)} product(s), only {diff} difference(s) "
                                  f"from a target amplicon (< min_snps={min_snps})")
            else:
                ev.reasons.append(f"background {bid}: {len(products)} product(s)")
        if i % 25 == 0 or i == len(background_files):
            log.info("  background %d/%d screened, %d pair(s) still passing", i, len(background_files),
                     sum(ev.passed for ev in alive))
    log.info("Background screen took %.0f s", time.time() - t_start)


def refine_loci(
    loci: Sequence[str],
    aln_by_locus: Mapping[str, Mapping[str, str]],
    genomes_by_id: Mapping[str, Genome],
    work: Path,
    threads: int,
) -> Dict[str, Tuple[Dict[str, str], LocusScore]]:
    """Re-align the ungapped sequences of several loci with mafft, `threads` loci at a time.

    Only loci that reach primer design are refined. A locus whose mafft run fails keeps
    its minimap2-based alignment."""
    (work / "unaligned").mkdir(parents=True, exist_ok=True)
    inputs = []
    for locus in loci:
        unaligned = work / "unaligned" / f"{locus}.fasta"
        write_fasta(unaligned, ((gid, s.replace("-", "")) for gid, s in aln_by_locus[locus].items()))
        inputs.append(unaligned)
    out: Dict[str, Tuple[Dict[str, str], LocusScore]] = {}
    for locus, (_, out_path) in zip(loci, align_loci(inputs, work, threads=threads)):
        if out_path is None:
            log.warning("mafft refinement failed for %s; using the minimap2-based alignment", locus)
            aln = dict(aln_by_locus[locus])
        else:
            aln = read_fasta_dict(out_path)
        out[locus] = (aln, score_alignment(locus, aln, genomes_by_id))
    return out


def _windows_overlap(a: Window, b: Window) -> bool:
    return a.contig == b.contig and a.start < b.end and b.start < a.end


def _reference_amplicon(ev: PairEvaluation, reference_id: str) -> Optional[Tuple[str, int, int]]:
    res = ev.per_genome.get(reference_id)
    if res is None or not res.on_target:
        return None
    pr = res.on_target[0]
    return pr.contig, pr.start, pr.end


def _regions_overlap(a: Optional[Tuple[str, int, int]], b: Optional[Tuple[str, int, int]]) -> bool:
    if a is None or b is None:
        return False
    return a[0] == b[0] and a[1] < b[2] and b[1] < a[2]


def _selection_rank_key(item: Tuple[LocusScore, CandidateLocus, PairEvaluation]) -> Tuple[bool, float, float]:
    """Best first: no background product, then the most informative locus, then the lowest primer3 penalty."""
    sc, _, ev = item
    return bool(ev.background_products), -sc.combined_score, ev.pair.penalty


def _rank_pair_table(df: pd.DataFrame) -> pd.DataFrame:
    """Selected pairs by Rank, then the other passing pairs, then failed ones; within those
    by locus score (high first) and primer3 penalty (low first)."""
    if df.empty:
        return df
    key = df.assign(_unselected=df["Rank"].isna(), _failed=~df["Passed_Checks"].astype(bool),
                    _bg=df["Background_Products"].fillna("").astype(str) != "")
    key = key.sort_values(["_unselected", "Rank", "_failed", "_bg", "Locus_Combined_Score", "Primer3_Penalty"],
                          ascending=[True, True, True, True, False, True], kind="stable")
    return df.loc[key.index].reset_index(drop=True)


def _pair_record(ev: PairEvaluation, score: LocusScore, rank: Optional[int], cfg: Config) -> Dict[str, object]:
    p = ev.pair
    rec: Dict[str, object] = {
        "Rank": rank,
        "Locus": p.locus,
        "Pair_Index": p.pair_index,
        "Selected": rank is not None,
        "Passed_Checks": ev.passed,
        "Failure_Reasons": " | ".join(ev.reasons),
        "LEFT_PRIMER": p.left_seq,
        "RIGHT_PRIMER": p.right_seq,
        "LEFT_TM": round(p.left_tm, 2),
        "RIGHT_TM": round(p.right_tm, 2),
        "Recommended_Annealing_Temp": round(cfg.recommended_annealing_temp(p.left_tm, p.right_tm), 1),
        "LEFT_GC": round(p.left_gc, 1),
        "RIGHT_GC": round(p.right_gc, 1),
        "Template_Product_Size": p.template_product_size,
        "Primer3_Penalty": round(p.penalty, 3),
        "Amplicon_Aln_Start": p.amplicon_start,
        "Amplicon_Aln_End": p.amplicon_end,
        "Locus_Combined_Score": round(score.combined_score, 5),
        "Informative_Sites_In_Amplicon_Total": sum(p.snps_in_amplicon.values()),
    }
    for key, n in sorted(p.snps_in_amplicon.items()):
        rec[f"Sites_In_Amplicon_{key}"] = n
    sizes, weakest = [], None
    for gid, res in ev.per_genome.items():
        sizes.append(f"{gid}:{','.join(str(pr.size) for pr in res.on_target) or '-'}")
        for pr in res.on_target:
            for tm in (pr.fwd_tm, pr.rev_tm):
                if tm is not None:
                    weakest = tm if weakest is None else min(weakest, tm)
    rec["On_Target_Product_Sizes"] = " ".join(sizes)
    rec["Background_Products"] = "; ".join(ev.background_products)
    rec["Min_On_Target_Site_Tm"] = round(weakest, 1) if weakest is not None else ""
    return rec


def amplicon_screen(
    selected: Sequence[Tuple[LocusScore, CandidateLocus, PairEvaluation]],
    genomes: Sequence[Genome],
    genome_seqs: Mapping[str, Mapping[str, str]],
    background_files: Sequence[Path],
    primer_dir: Path,
    cfg: Config,
) -> Optional[pd.DataFrame]:
    """Similarity of the selected reference amplicons to background genomes and (optionally) NCBI nt."""
    from . import blast

    asp = cfg.amplicon_screen
    if not (asp.background_identity and background_files) and not asp.remote_blast:
        return None
    out_dir = primer_dir / "amplicon_screen"
    out_dir.mkdir(exist_ok=True)
    reference_id = choose_reference(genomes).genome_id
    records, tolerated, amp_len = [], {}, {}
    for sc, cand, ev in selected:
        res = ev.per_genome.get(reference_id)
        if res is None or not res.on_target:
            continue
        pr = res.on_target[0]
        seq = product_sequence(genome_seqs[reference_id], pr)
        records.append((sc.locus, seq))
        amp_len[sc.locus] = len(seq)
        ext = [v for k, v in ev.pair.snps_in_amplicon.items() if k.startswith("External_")]
        tolerated[sc.locus] = blast.tolerated_identity(len(seq), min(ext) if ext else None, asp.fallback_tolerated_identity)
    if not records:
        return None
    amplicons_fa = out_dir / "reference_amplicons.fasta"
    write_fasta(amplicons_fa, records)
    target_species = sorted({g.species for g in genomes if g.is_target})
    nontarget_species = sorted({g.species for g in genomes if not g.is_target})
    summary = pd.DataFrame({"Locus": [r[0] for r in records], "Amplicon_Length": [amp_len[r[0]] for r in records],
                            "Tolerated_Identity": [round(tolerated[r[0]], 2) for r in records]})

    if asp.background_identity and background_files:
        bg = blast.background_identity(amplicons_fa, background_files, out_dir / "tmp", cfg.max_cores, tolerated)
        shutil.rmtree(out_dir / "tmp", ignore_errors=True)
        bg.to_csv(out_dir / "background_identity.csv", index=False)
        best = bg.sort_values("Pct_Identity", ascending=False).drop_duplicates("Locus")
        summary = summary.merge(best[["Locus", "Background_Genome", "Pct_Identity"]].rename(
            columns={"Background_Genome": "Closest_Background_Genome", "Pct_Identity": "Background_Max_Identity"}), how="left")
        summary["Background_Risk_Genomes"] = summary["Locus"].map(
            lambda l: ";".join(bg[(bg.Locus == l) & bg.Risk]["Background_Genome"]))
        for _, r in summary[summary["Background_Risk_Genomes"].astype(bool)].iterrows():
            log.warning("%s: background genome(s) %s are at least as similar to the amplicon as the nearest non-target",
                        r["Locus"], r["Background_Risk_Genomes"])

    if asp.remote_blast:
        try:
            xml = blast.remote_blast(records, asp.blast_db, asp.megablast, asp.hitlist_size, asp.entrez_query)
            (out_dir / "blast_results.xml").write_text(xml)
            from io import StringIO
            hits = blast.parse_blast_xml(StringIO(xml), asp.min_identity_report)
            table = blast.summarize_blast(hits, target_species, nontarget_species, tolerated)
            table.to_csv(out_dir / "blast_hits.csv", index=False)
            summary["BLAST_N_Hits"] = summary["Locus"].map(lambda l: int((table.Locus == l).sum()))
            summary["BLAST_Unknown_Taxa"] = summary["Locus"].map(
                lambda l: ";".join(sorted(set(table[(table.Locus == l) & (table.Status == blast.STATUS_UNKNOWN)]["Organism"]))))
            summary["BLAST_Risk_Taxa"] = summary["Locus"].map(
                lambda l: ";".join(sorted(set(table[(table.Locus == l) & table.Risk]["Organism"]))))
            for _, r in summary[summary["BLAST_Risk_Taxa"].astype(bool)].iterrows():
                log.warning("%s: NCBI nt organisms %s are at least as similar to the amplicon as the nearest non-target; "
                            "add them as non-targets or background", r["Locus"], r["BLAST_Risk_Taxa"])
        except Exception as e:  # network / NCBI problems must not kill the run
            log.error("Remote BLAST failed: %s", e)
            summary["BLAST_Error"] = str(e)

    summary.to_csv(out_dir / "amplicon_screen_summary.csv", index=False)
    return summary


def run_group(group: str, genomes: Sequence[Genome], cfg: Config) -> Dict[str, object]:
    t0 = time.time()
    out = cfg.output_dir / group
    out.mkdir(parents=True, exist_ok=True)
    genomes_by_id = {g.genome_id: g for g in genomes}
    log.info("Group %s: %d target genome(s), %d non-target genome(s)",
             group, sum(g.is_target for g in genomes), sum(not g.is_target for g in genomes))

    reference = choose_reference(genomes)
    candidates, genome_seqs = discover_loci(genomes, reference, out / "locus_discovery", cfg.locus_discovery, cfg.max_cores)
    log.info("Locus discovery took %.0f s", time.time() - t0)
    if not candidates:
        raise RuntimeError("No single-copy loci shared by the required genomes were found")

    # ---- alignments from the minimap2 cs tags; no MSA program needed -------------
    aligned_dir = out / "loci" / "aligned"
    aligned_dir.mkdir(parents=True, exist_ok=True)
    ref_seqs = genome_seqs[reference.genome_id]
    scores: List[LocusScore] = []
    cand_by_locus: Dict[str, CandidateLocus] = {}
    aln_by_locus: Dict[str, Dict[str, str]] = {}
    locus_rows = []
    for cand in candidates:
        locus = locus_id(cand.window)
        w = cand.window
        for gid, h in cand.hits.items():
            locus_rows.append({"Locus": locus, "Reference_Contig": w.contig, "Reference_Start": w.start,
                               "Reference_End": w.end, "Genome": gid, "Contig": h.contig, "Start": h.start,
                               "End": h.end, "Strand": h.strand, "Divergence": h.divergence})
        window_seq = ref_seqs[w.contig][w.start : w.end]
        try:
            aln = build_alignment(window_seq, {gid: (h.qstart, h.qend, h.strand, h.cs) for gid, h in cand.hits.items()})
            sc = score_alignment(locus, aln, genomes_by_id)
        except ValueError as e:
            log.warning("Skipping %s: %s", locus, e)
            continue
        write_fasta(aligned_dir / f"{locus}.fasta", aln.items())
        scores.append(sc)
        cand_by_locus[locus] = cand
        aln_by_locus[locus] = aln
    pd.DataFrame(locus_rows).to_csv(out / "loci" / "locus_table.csv", index=False)

    info_dir = out / "informative_loci"
    info_dir.mkdir(exist_ok=True)
    summary_csv = info_dir / "snp_summary.csv"
    df = pd.DataFrame([s.to_record() for s in scores])
    if not df.empty:
        df.sort_values("Combined_Informativeness_Score", ascending=False, inplace=True)
    df.to_csv(summary_csv, index=False)
    for sc in scores:
        if sc.combined_score >= cfg.snp_avg_prop_threshold:
            shutil.copy(aligned_dir / f"{sc.locus}.fasta", info_dir / f"{sc.locus}.fasta")
    log.info("Scored %d loci; %d pass snp_avg_prop_threshold=%.3g",
             len(scores), int((df["Combined_Informativeness_Score"] >= cfg.snp_avg_prop_threshold).sum()) if not df.empty else 0,
             cfg.snp_avg_prop_threshold)

    if cfg.make_plots and not df.empty:
        from .plots import plot_informativeness_heatmap
        plot_informativeness_heatmap(summary_csv, info_dir / "heatmap_informativeness.png")

    # ---- primer design and in silico screening ----------------------------------
    primer_dir = out / "primers"
    primer_dir.mkdir(exist_ok=True)
    cons_dir = info_dir / "consensus_sequences"
    cons_dir.mkdir(exist_ok=True)
    amplicon_dir = primer_dir / "amplicons"
    amplicon_dir.mkdir(exist_ok=True)

    sel = cfg.primer_selection
    binding = binding_params_from_config(cfg)
    if binding.min_tm is not None:
        log.info("Binding model: <= %d mismatch(es) in the 3'-terminal %d nt, <= %d overall, duplex Tm >= %.1f C "
                 "(annealing %.1f C minus %.1f C margin; primer %.0f nM)", binding.three_prime_max_mismatches,
                 binding.three_prime_region, binding.max_mismatches, binding.min_tm, cfg.annealing_temp,
                 sel.binding_tm_margin, float(cfg.primer3_params.get("PRIMER_DNA_CONC", 250.0)))
        if binding.near_miss_tm is not None:
            log.info("Near misses: products whose sites reach only %.1f-%.1f C also reject a pair",
                     binding.near_miss_tm, binding.min_tm)
    else:
        log.info("Binding model: mismatch counts only, <= %d in the 3'-terminal %d nt and <= %d overall",
                 binding.three_prime_max_mismatches, binding.three_prime_region, binding.max_mismatches)
    background_files = list_background(sel.background_dir)
    if background_files:
        log.info("%d background genome(s) in %s", len(background_files), sel.background_dir)
    else:
        log.warning("primer_selection.background_dir is not set: primers are screened only against the %d input "
                    "genomes. Off-target products from other organisms in your samples cannot be detected. "
                    "Add genomes of host, contaminants and expected community members as background.", len(genomes))
    indexes = {gid: GenomeIndex(seqs) for gid, seqs in genome_seqs.items()}
    if _fork_pool_possible(cfg.max_cores):
        # Build the k-mer lookups before any worker is forked, so all workers share them.
        # numpy releases the GIL while sorting, so threads build them in parallel.
        with ThreadPoolExecutor(max_workers=cfg.max_cores) as tp:
            list(tp.map(lambda ix: ix.warm(_seed_k(binding)), indexes.values()))
        _SHARED.update(genomes=genomes, indexes=indexes, cands=cand_by_locus, cfg=cfg, binding=binding)

    p3 = dict(cfg.primer3_params)
    sd = cfg.snp_primer_design
    max_target_span = cfg.max_product_size - 2 * int(p3["PRIMER_MAX_SIZE"]) - 20
    ranked = sorted(scores, key=lambda s: -s.combined_score)
    locus_records: Dict[str, Dict[str, object]] = {}
    designed: List[Tuple[LocusScore, CandidateLocus, List[PairEvaluation]]] = []
    passing_windows: List[Window] = []
    selected: List[Tuple[LocusScore, CandidateLocus, PairEvaluation]] = []
    refined: Dict[str, Tuple[Dict[str, str], LocusScore]] = {}
    n_try = min(sd.max_loci_to_try, len(ranked))
    next_rank = 0

    reserve: List[Tuple[LocusScore, CandidateLocus, PairEvaluation]] = []
    from . import multiplex
    # Cross-pair products of the selected pairs; near misses count here too.
    pool_binding = binding if binding.near_miss_tm is None else replace(binding, min_tm=binding.near_miss_tm,
                                                                        near_miss_tm=None)

    def select(sc: LocusScore, cand: CandidateLocus, chosen: PairEvaluation, note: str) -> None:
        if len(selected) >= sd.snp_top_n:
            locus_records[sc.locus]["Outcome"] += "; not selected, snp_top_n already reached"
            return
        # Overlapping windows often hold the same amplicon; compare amplicons on the reference.
        amp = _reference_amplicon(chosen, reference.genome_id)
        overlapping = [s.locus for s, _, e in selected if _regions_overlap(amp, _reference_amplicon(e, reference.genome_id))]
        if overlapping:
            locus_records[sc.locus]["Outcome"] = f"skipped: amplicon overlaps selected locus {overlapping[0]}"
            return
        if sel.poolable and selected:
            # the primers of this pair must not form products (or near misses) with any selected pair
            cross = multiplex.cross_pair_products(
                [(s.locus, e.pair.left_seq, e.pair.right_seq) for s, _, e in selected]
                + [(sc.locus, chosen.pair.left_seq, chosen.pair.right_seq)],
                indexes, pool_binding, sel.max_offtarget_product_size)
            partners = sorted({b if a == sc.locus else a for a, b in multiplex.conflicting_loci(cross) if sc.locus in (a, b)})
            if partners:
                locus_records[sc.locus]["Outcome"] = (f"skipped: primers form cross-pair products with selected locus "
                                                      f"{partners[0]} (poolable: true)")
                return
        locus_records[sc.locus]["Outcome"] = f"selected pair {chosen.pair.pair_index}{note}"
        selected.append((sc, cand, chosen))
        _write_selected(sc, cand, chosen, genome_seqs, primer_dir, amplicon_dir)

    # Rounds of: (A) design and screen loci against the group's own genomes, (B) one
    # background pass for the new survivors, (C) select. A round designs a few more loci
    # than still missing when a background screen may reject some; further rounds run
    # until snp_top_n loci are selected or max_loci_to_try loci have been tried.
    while len(selected) < sd.snp_top_n and next_rank < n_try:
        missing = sd.snp_top_n - len(selected)
        wanted = missing * (2 if background_files else 1)
        new_passing = 0
        round_designed: List[Tuple[LocusScore, CandidateLocus, List[PairEvaluation]]] = []

        # Phase A: primer3 for a batch of loci, in silico PCR of all their pairs (one pair
        # per worker process), then bookkeeping in rank order. The walk stops exactly where
        # a one-locus-at-a-time loop would stop, so results do not depend on max_cores; a
        # parallel batch only wastes the in silico PCR of the loci after the stop point.
        t_round, rank_before = time.time(), next_rank
        pool = _fork_pool(cfg.max_cores)
        batch_loci = 1 if pool is None else max(1, -(-cfg.max_cores // max(1, int(p3.get("PRIMER_NUM_RETURN", 5)))))
        try:
            while next_rank < n_try and new_passing < wanted:
                designs: List[_LocusDesign] = []
                for r in range(next_rank, min(n_try, next_rank + batch_loci)):
                    sc = ranked[r]
                    aln = aln_by_locus[sc.locus]
                    if cfg.locus_discovery.refine_with_mafft:
                        if sc.locus not in refined:
                            # refine the next max_cores loci at once; mafft runs are independent
                            todo = [s.locus for s in ranked[r : min(n_try, r + cfg.max_cores)] if s.locus not in refined]
                            refined.update(refine_loci(todo, aln_by_locus, genomes_by_id, out / "loci" / "mafft", cfg.max_cores))
                        aln, sc = refined[sc.locus]
                    designs.append(_design_pairs(sc, aln, cand_by_locus[sc.locus], genomes_by_id, cfg, p3, max_target_span))
                jobs = [(pair, d.sc.locus) for d in designs for pair in d.pairs]
                if pool is None:
                    evs_flat = [evaluate_pair(pair, genomes, indexes, cand_by_locus[locus], cfg, binding) for pair, locus in jobs]
                else:
                    evs_flat = list(pool.map(_evaluate_job, jobs))
                pos = 0
                for d in designs:
                    if new_passing >= wanted:
                        break
                    next_rank += 1
                    evs = evs_flat[pos : pos + len(d.pairs)]
                    pos += len(d.pairs)
                    ev_list = _record_design(d, evs, cons_dir, locus_records)
                    if ev_list is None:
                        continue
                    round_designed.append((d.sc, d.cand, ev_list))
                    if any(ev.passed for ev in ev_list) and not any(_windows_overlap(d.cand.window, w) for w in passing_windows):
                        passing_windows.append(d.cand.window)
                        new_passing += 1
        finally:
            if pool is not None:
                pool.shutdown()
        designed.extend(round_designed)
        log.info("Design round: %d locus/loci tried, %d pair(s) tested in silico, %d new passing region(s), %.0f s",
                 next_rank - rank_before, sum(len(e) for _, _, e in round_designed), new_passing, time.time() - t_round)

        # Phase B
        survivors = [ev for _, _, evs in round_designed for ev in evs if ev.passed]
        if background_files and survivors:
            for ev in survivors:
                ev.target_amplicons = sorted({product_sequence(genome_seqs[g.genome_id], pr).upper()
                                              for g in genomes if g.is_target
                                              for pr in ev.per_genome[g.genome_id].on_target})
            screen_background(survivors, background_files, binding, sel.max_offtarget_product_size, cfg.max_cores,
                              sel.background_products, sd.min_snps)

        # Phase C: the best surviving pair per locus, one locus per region, in rank order.
        # Pairs without any background product are preferred; a locus whose passing pairs
        # all amplify in the background (distinguishably) is kept in reserve.
        for sc, cand, evs in round_designed:
            passing = [ev for ev in evs if ev.passed]
            if not passing:
                if any(r.startswith("background") for ev in evs for r in ev.reasons):
                    locus_records[sc.locus]["Outcome"] = "all passing pair(s) rejected by the background screen"
                continue
            clean = next((ev for ev in passing if not ev.background_products), None)
            if clean is None:
                reserve.append((sc, cand, passing[0]))
                locus_records[sc.locus]["Outcome"] = ("passing pair(s) amplify in the background, distinguishable "
                                                      "by sequence; kept in reserve")
                continue
            select(sc, cand, clean, "")

    # Fallback: too few loci without background products. Fill the remaining slots with
    # loci whose background products differ from the target amplicons by >= min_snps.
    if len(selected) < sd.snp_top_n and reserve:
        log.warning("Group %s: only %d locus/loci without background products; adding loci whose background "
                    "products differ from the target amplicons by >= min_snps=%d (see Background_Products)",
                    group, len(selected), sd.min_snps)
        for sc, cand, ev in reserve:
            select(sc, cand, ev, " (fallback: amplifies in the background, distinguishable by sequence)")

    _SHARED.clear()  # release the indexes before the next group

    # Every output lists the selected pairs best first (Rank 1 = best).
    selected.sort(key=_selection_rank_key)
    rank_of = {id(ev): i for i, (_, _, ev) in enumerate(selected, 1)}
    for i, (sc, _, ev) in enumerate(selected, 1):
        log.info("Rank %d: %s pair %d (locus score %.4f, primer3 penalty %.3f%s)", i, sc.locus, ev.pair.pair_index,
                 sc.combined_score, ev.pair.penalty, ", amplifies in the background" if ev.background_products else "")
    write_fasta(primer_dir / "selected_primers.fasta",
                [rec for i, (sc, _, ev) in enumerate(selected, 1)
                 for rec in ((f"rank{i}_{sc.locus}_LEFT", ev.pair.left_seq), (f"rank{i}_{sc.locus}_RIGHT", ev.pair.right_seq))])

    # Loci that amplify in the background: show how their amplicons differ from the target.
    if cfg.make_plots:
        for sc, cand, ev in selected:
            if ev.background_sequences:
                _plot_alignment(sc, ev, genomes, genome_seqs, reference.genome_id, primer_dir / "amplicon_alignments",
                                sd.min_snps)

    # Reference panel and diagnostic sites for PRISCREEN (read classification).
    if selected:
        _write_priscreen_exports(selected, genomes, genome_seqs, reference.genome_id, primer_dir)

    # Pooling the selected pairs in one reaction: products formed by primers of two pairs
    # (none with poolable: true).
    cross = multiplex.cross_pair_products([(sc.locus, ev.pair.left_seq, ev.pair.right_seq) for sc, _, ev in selected],
                                          indexes, pool_binding, sel.max_offtarget_product_size)
    cross.to_csv(primer_dir / "multiplex_check.csv", index=False)
    for a, b in multiplex.conflicting_loci(cross):
        log.warning("Selected loci %s and %s form cross-pair products; do not pool them in one reaction "
                    "(see primers/multiplex_check.csv)", a, b)

    if len(selected) < sd.snp_top_n:
        log.warning("Group %s: only %d of snp_top_n=%d loci selected after trying %d locus/loci "
                    "(max_loci_to_try=%d, %d scored)", group, len(selected), sd.snp_top_n, next_rank,
                    sd.max_loci_to_try, len(ranked))

    screen_summary = amplicon_screen(selected, genomes, genome_seqs, background_files, primer_dir, cfg) if selected else None

    pair_records = [_pair_record(ev, sc, rank_of.get(id(ev)), cfg) for sc, _, evs in designed for ev in evs]
    all_pairs = _rank_pair_table(pd.DataFrame(pair_records))
    if "Rank" in all_pairs:
        all_pairs["Rank"] = all_pairs["Rank"].astype("Int64")  # blank, not NaN, for unselected pairs
    all_pairs.to_csv(primer_dir / "primer_design_summary.csv", index=False)
    pd.DataFrame(list(locus_records.values())).to_csv(primer_dir / "locus_design_log.csv", index=False)
    # keep the header even when nothing was selected, so the file is still a valid table
    sel_df = all_pairs[all_pairs["Selected"]] if "Selected" in all_pairs else all_pairs
    sel_df.to_csv(primer_dir / "selected_primers.csv", index=False)

    if cfg.make_plots and selected:
        from .plots import plot_snp_positions
        plot_snp_positions([sc for sc, _, _ in selected], info_dir / "snp_density_plots")

    elapsed = time.time() - t0
    log.info("Group %s: %d primer set(s) selected in %.0f s", group, len(selected), elapsed)
    return {
        "Group": group,
        "Reference": reference.genome_id,
        "N_Genomes": len(genomes),
        "N_Background_Genomes": len(background_files),
        "N_Candidate_Loci": len(candidates),
        "N_Scored_Loci": len(scores),
        "N_Loci_Tried": next_rank,
        "N_Loci_Designed": len(designed),
        "N_Selected_Primer_Sets": len(selected),
        "Seconds": round(elapsed),
    }


@dataclass
class _LocusDesign:
    """Primer3 result for one locus, before in silico PCR and before anything is written."""
    sc: LocusScore
    cand: CandidateLocus
    record: Dict[str, object]
    consensus: str
    template: str
    pairs: List[PrimerPair] = field(default_factory=list)


def _design_pairs(
    sc: LocusScore,
    aln: Mapping[str, str],
    cand: CandidateLocus,
    genomes_by_id: Mapping[str, Genome],
    cfg: Config,
    p3: Mapping[str, object],
    max_target_span: int,
) -> _LocusDesign:
    """Primer3 design for one locus. No pairs means none could be designed; the reason is in
    `record["Outcome"]`."""
    sel = cfg.primer_selection
    sd = cfg.snp_primer_design
    cons_ids = list(aln) if sel.conservation == "all" else [gid for gid in aln if genomes_by_id[gid].is_target]
    template = masked_template([aln[g] for g in cons_ids])
    lrec: Dict[str, object] = {"Locus": sc.locus, "Combined_Score": sc.combined_score,
                               "Conserved_Fraction": 1 - template.count("N") / len(template)}
    d = _LocusDesign(sc, cand, lrec, iupac_consensus(list(aln.values())), template)
    target = choose_target_region(sc.positions_by_key(), sd.snp_window_size, sd.min_snps, max_target_span, len(template))
    if target is None:
        lrec["Outcome"] = f"no window of {sd.snp_window_size} bp with >= {sd.min_snps} informative sites"
        return d
    lrec["Target_Start"], lrec["Target_Length"] = target
    try:
        pairs, explain = design_primers(sc.locus, template, target, p3)
    except (OSError, RuntimeError) as e:
        lrec["Outcome"] = f"primer3 error: {e}"
        return d
    if not pairs:
        lrec["Outcome"] = f"primer3 found no pair ({explain})"
        return d
    for pair in pairs:
        pair.snps_in_amplicon = sc.count_in_region(pair.amplicon_start, pair.amplicon_end)
    d.pairs = list(pairs)
    return d


def _record_design(
    d: _LocusDesign,
    evs: Sequence[PairEvaluation],
    cons_dir: Path,
    locus_records: Dict[str, Dict[str, object]],
) -> Optional[List[PairEvaluation]]:
    """Write the consensus files and the log record of a locus the design walk reached.

    Returns the pair evaluations, or None when no pair could be designed."""
    locus = d.sc.locus
    write_fasta(cons_dir / f"{locus}.fasta", [(f"{locus} IUPAC consensus (gaps ignored)", d.consensus)])
    write_fasta(cons_dir / f"{locus}.masked_template.fasta", [(f"{locus} N-masked design template", d.template)])
    locus_records[locus] = d.record
    if not d.pairs:
        return None
    n_pass = sum(ev.passed for ev in evs)
    if n_pass:
        d.record["Outcome"] = f"{n_pass} of {len(d.pairs)} pair(s) passed on the input genomes"
    else:
        d.record["Outcome"] = f"{len(d.pairs)} pair(s) designed, none passed in silico checks"
    return list(evs)


def _plot_alignment(sc: LocusScore, ev: PairEvaluation, genomes: Sequence[Genome],
                    genome_seqs: Mapping[str, Mapping[str, str]], reference_id: str, out_dir: Path,
                    min_snps: int) -> None:
    """Alignment figure (and difference table) of the target, non-target and background
    amplicons of a selected pair, against the reference amplicon."""
    from .plots import plot_amplicon_alignment

    ref_res = ev.per_genome.get(reference_id)
    if ref_res is None or not ref_res.on_target:
        return
    rows = []
    for g in genomes:
        res = ev.per_genome.get(g.genome_id)
        for pr in (res.on_target if res else []):
            rows.append((g.genome_id, "target" if g.is_target else "non-target", product_sequence(genome_seqs[g.genome_id], pr)))
    rows += [(bid, "background", seq) for bid, seq in sorted(ev.background_sequences.items())]
    out_dir.mkdir(parents=True, exist_ok=True)
    note = (f"background genomes amplify; each background product differs from every target amplicon "
            f"by >= {min_snps} substitutions + indel events")
    plot_amplicon_alignment(sc.locus, (reference_id, product_sequence(genome_seqs[reference_id], ref_res.on_target[0])),
                            rows, out_dir / f"{sc.locus}.png", out_dir / f"{sc.locus}_differences.tsv", note)


def _panel_rows(cand: CandidateLocus, ev: PairEvaluation, genomes: Sequence[Genome],
                genome_seqs: Mapping[str, Mapping[str, str]], reference_id: str):
    """(reference amplicon, [(genome_id, species, role, amplicon)]) of a selected pair, all
    in the orientation of the reference locus; None without a reference amplicon."""
    from . import panel

    ref_res = ev.per_genome.get(reference_id)
    if ref_res is None or not ref_res.on_target:
        return None
    rows = []
    for g in genomes:
        res = ev.per_genome.get(g.genome_id)
        hit = cand.hits.get(g.genome_id)
        for pr in (res.on_target if res else []):
            seq = product_sequence(genome_seqs[g.genome_id], pr).upper()
            if hit is not None and hit.strand == "-":
                seq = revcomp(seq)
            rows.append((g.genome_id, g.species, "target" if g.is_target else "non-target", seq))
    reference = next(seq for gid, _, _, seq in rows if gid == reference_id)
    for bid, seq in sorted(ev.background_sequences.items()):
        rows.append((bid, species_from_stem(bid), "background", panel.orient(reference, seq)))
    return reference, rows


def _write_priscreen_exports(selected, genomes: Sequence[Genome], genome_seqs: Mapping[str, Mapping[str, str]],
                             reference_id: str, primer_dir: Path) -> None:
    """primers/priscreen_panel.fasta and primers/diagnostic_sites/<locus>.tsv."""
    from . import panel

    loci = []
    for sc, cand, ev in selected:
        got = _panel_rows(cand, ev, genomes, genome_seqs, reference_id)
        if got is None:
            continue
        reference, rows = got
        loci.append((sc.locus, rows))
        n = panel.write_diagnostic_sites(primer_dir / "diagnostic_sites" / f"{sc.locus}.tsv", reference, rows)
        log.info("Locus %s: %d diagnostic site(s) written", sc.locus, n)
    panel.write_panel(primer_dir / "priscreen_panel.fasta", loci)


def _write_selected(
    sc: LocusScore,
    cand: CandidateLocus,
    chosen: PairEvaluation,
    genome_seqs: Mapping[str, Mapping[str, str]],
    primer_dir: Path,
    amplicon_dir: Path,
) -> None:
    """Primer FASTA and predicted amplicons of a selected pair."""
    p = chosen.pair
    write_fasta(primer_dir / f"{sc.locus}_primers.fasta",
                [(f"{sc.locus}_LEFT", p.left_seq), (f"{sc.locus}_RIGHT", p.right_seq)])
    # Amplicons are written in the orientation of the reference locus, so that
    # they can be aligned or used as a reference database directly.
    amp_records = []
    for gid, res in chosen.per_genome.items():
        hit = cand.hits.get(gid)
        for pr in res.on_target:
            seq = product_sequence(genome_seqs[gid], pr)
            if hit is not None and hit.strand == "-":
                seq = revcomp(seq)
            amp_records.append((f"{gid} {pr.contig}:{pr.start}-{pr.end} strand={hit.strand if hit else '?'} size={pr.size}", seq))
    write_fasta(amplicon_dir / f"{sc.locus}_amplicons.fasta", amp_records)


def run_pipeline(cfg: Config) -> pd.DataFrame:
    log.info("PRISTINE %s | input: %s | output: %s | threads: %d", __version__, cfg.input_dir, cfg.output_dir, cfg.max_cores)
    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    groups = discover_groups(cfg.input_dir)
    rows = []
    for group, genomes in groups.items():
        log.info("=== %s ===", group)
        try:
            rows.append(run_group(group, genomes, cfg))
        except Exception as e:  # keep going with the other groups, but be loud
            log.error("Group %s failed: %s\n%s", group, e, traceback.format_exc())
            rows.append({"Group": group, "Error": str(e)})
    from .panel import merge_panels
    merge_panels({g: cfg.output_dir / g / "primers" / "priscreen_panel.fasta" for g in groups},
                 cfg.output_dir / "priscreen_panel.fasta")
    df = pd.DataFrame(rows)
    df.to_csv(cfg.output_dir / "run_summary.csv", index=False)
    return df
