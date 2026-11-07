"""Orchestration of the PRISTINE pipeline for every target group."""

from __future__ import annotations

import logging
import shutil
import time
import traceback
from dataclasses import dataclass, field
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
    Thermo,
    amplify,
    classify_products,
    product_sequence,
)
from .loci import CandidateLocus, Window, discover_loci, locus_id
from .manifest import Genome, choose_reference, discover_groups, is_fasta
from .primers import PrimerPair, choose_target_region, design_primers
from .pseudoalign import build_alignment
from .scoring import LocusScore, score_alignment
from .seqio import read_fasta_dict, revcomp, write_fasta

log = logging.getLogger(__name__)


@dataclass
class PairEvaluation:
    pair: PrimerPair
    passed: bool
    reasons: List[str] = field(default_factory=list)
    per_genome: Dict[str, GenomePcrResult] = field(default_factory=dict)


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
            dna_conc=float(p3.get("PRIMER_DNA_CONC", 50.0)),
        )
    return BindingParams(max_mismatches=sel.max_mismatches, three_prime_region=sel.three_prime_region,
                         three_prime_max_mismatches=sel.three_prime_max_mismatches, seed_length=sel.seed_length,
                         min_tm=cfg.min_site_tm, thermo=thermo)


def evaluate_pair(
    pair: PrimerPair,
    genomes: Sequence[Genome],
    indexes: Mapping[str, GenomeIndex],
    locus: CandidateLocus,
    cfg: Config,
    binding: BindingParams,
) -> PairEvaluation:
    """In silico PCR of one pair against the group's genomes plus the informative-site criterion."""
    sel = cfg.primer_selection
    ev = PairEvaluation(pair=pair, passed=True)
    size_lo = int(cfg.min_product_size * 0.9)
    size_hi = int(cfg.max_product_size * 1.1)

    for g in genomes:
        hit = locus.hits.get(g.genome_id)
        region = (hit.contig, hit.start, hit.end) if hit else None
        products = amplify(indexes[g.genome_id], pair.left_seq, pair.right_seq, sel.max_offtarget_product_size, binding)
        res = classify_products(g.genome_id, products, region)
        ev.per_genome[g.genome_id] = res
        if res.off_target:
            ev.passed = False
            ev.reasons.append(f"{g.genome_id}: {len(res.off_target)} off-target product(s)")
        n_on = len(res.on_target)
        if (g.is_target or sel.conservation == "all") and n_on == 0:
            ev.passed = False
            ev.reasons.append(f"{g.genome_id}: no product at locus")
        if n_on > 1:
            ev.passed = False
            ev.reasons.append(f"{g.genome_id}: {n_on} products at locus")
        for pr in res.on_target:
            if not (size_lo <= pr.size <= size_hi):
                ev.passed = False
                ev.reasons.append(f"{g.genome_id}: product size {pr.size} outside {size_lo}-{size_hi}")

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
    return ev


def screen_background(
    evaluations: Sequence[PairEvaluation],
    background_files: Sequence[Path],
    binding: BindingParams,
    max_product_size: int,
) -> None:
    """One streaming pass over the background genomes; every product there fails the pair."""
    alive = [ev for ev in evaluations if ev.passed]
    log.info("Screening %d candidate pair(s) against %d background genome(s)", len(alive), len(background_files))
    for i, path in enumerate(background_files, 1):
        alive = [ev for ev in alive if ev.passed]
        if not alive:
            break
        index = GenomeIndex(read_fasta_dict(path))
        bid = path.stem
        for ev in alive:
            products = amplify(index, ev.pair.left_seq, ev.pair.right_seq, max_product_size, binding)
            if products:
                ev.passed = False
                ev.reasons.append(f"background {bid}: {len(products)} product(s)")
                ev.per_genome[bid] = classify_products(bid, products, None)
        if i % 25 == 0 or i == len(background_files):
            log.info("  background %d/%d screened, %d pair(s) still passing", i, len(background_files),
                     sum(ev.passed for ev in alive))


def refine_locus(locus: str, aln: Mapping[str, str], genomes_by_id: Mapping[str, Genome], work: Path):
    """Re-align the ungapped locus sequences with mafft (only for loci that reach primer design)."""
    (work / "unaligned").mkdir(parents=True, exist_ok=True)
    unaligned = work / "unaligned" / f"{locus}.fasta"
    write_fasta(unaligned, ((gid, s.replace("-", "")) for gid, s in aln.items()))
    (_, out_path), = align_loci([unaligned], work, threads=1)
    if out_path is None:
        log.warning("mafft refinement failed for %s; using the minimap2-based alignment", locus)
        return dict(aln), score_alignment(locus, aln, genomes_by_id)
    refined = read_fasta_dict(out_path)
    return refined, score_alignment(locus, refined, genomes_by_id)


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


def _pair_record(ev: PairEvaluation, score: LocusScore, selected: bool) -> Dict[str, object]:
    p = ev.pair
    rec: Dict[str, object] = {
        "Locus": p.locus,
        "Pair_Index": p.pair_index,
        "Selected": selected,
        "Passed_Checks": ev.passed,
        "Failure_Reasons": " | ".join(ev.reasons),
        "LEFT_PRIMER": p.left_seq,
        "RIGHT_PRIMER": p.right_seq,
        "LEFT_TM": round(p.left_tm, 2),
        "RIGHT_TM": round(p.right_tm, 2),
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
                 "(annealing %.1f C minus %.1f C margin)", binding.three_prime_max_mismatches, binding.three_prime_region,
                 binding.max_mismatches, binding.min_tm, cfg.annealing_temp, sel.binding_tm_margin)
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

    p3 = dict(cfg.primer3_params)
    sd = cfg.snp_primer_design
    max_target_span = cfg.max_product_size - 2 * int(p3["PRIMER_MAX_SIZE"]) - 20
    ranked = sorted(scores, key=lambda s: -s.combined_score)
    locus_records: Dict[str, Dict[str, object]] = {}
    designed: List[Tuple[LocusScore, CandidateLocus, List[PairEvaluation]]] = []
    passing_windows: List[Window] = []
    # Design a few more loci than requested when a background screen may still reject some.
    wanted = sd.snp_top_n * (2 if background_files else 1)

    # Phase A: design and screen against the group's own genomes.
    tried = 0
    for sc in ranked:
        if tried >= sd.max_loci_to_try or len(passing_windows) >= wanted:
            break
        tried += 1
        aln = aln_by_locus[sc.locus]
        cand = cand_by_locus[sc.locus]
        if cfg.locus_discovery.refine_with_mafft:
            aln, sc = refine_locus(sc.locus, aln, genomes_by_id, out / "loci" / "mafft")
        cons_ids = list(aln) if sel.conservation == "all" else [gid for gid in aln if genomes_by_id[gid].is_target]
        template = masked_template([aln[g] for g in cons_ids])
        write_fasta(cons_dir / f"{sc.locus}.fasta", [(f"{sc.locus} IUPAC consensus (gaps ignored)", iupac_consensus(list(aln.values())))])
        write_fasta(cons_dir / f"{sc.locus}.masked_template.fasta", [(f"{sc.locus} N-masked design template", template)])

        lrec: Dict[str, object] = {"Locus": sc.locus, "Combined_Score": sc.combined_score,
                                   "Conserved_Fraction": 1 - template.count("N") / len(template)}
        locus_records[sc.locus] = lrec
        target = choose_target_region(sc.positions_by_key(), sd.snp_window_size, sd.min_snps, max_target_span, len(template))
        if target is None:
            lrec["Outcome"] = f"no window of {sd.snp_window_size} bp with >= {sd.min_snps} informative sites"
            continue
        lrec["Target_Start"], lrec["Target_Length"] = target
        try:
            pairs, explain = design_primers(sc.locus, template, target, p3)
        except (OSError, RuntimeError) as e:
            lrec["Outcome"] = f"primer3 error: {e}"
            continue
        if not pairs:
            lrec["Outcome"] = f"primer3 found no pair ({explain})"
            continue
        evs = []
        for pair in pairs:
            pair.snps_in_amplicon = sc.count_in_region(pair.amplicon_start, pair.amplicon_end)
            evs.append(evaluate_pair(pair, genomes, indexes, cand, cfg, binding))
        designed.append((sc, cand, evs))
        n_pass = sum(ev.passed for ev in evs)
        if n_pass:
            lrec["Outcome"] = f"{n_pass} of {len(pairs)} pair(s) passed on the input genomes"
            if not any(_windows_overlap(cand.window, w) for w in passing_windows):
                passing_windows.append(cand.window)
        else:
            lrec["Outcome"] = f"{len(pairs)} pair(s) designed, none passed in silico checks"

    # Phase B: one pass over the background genomes for all surviving pairs.
    survivors = [ev for _, _, evs in designed for ev in evs if ev.passed]
    if background_files and survivors:
        screen_background(survivors, background_files, binding, sel.max_offtarget_product_size)

    # Phase C: pick the best surviving pair per locus, one locus per region.
    selected: List[Tuple[LocusScore, CandidateLocus, PairEvaluation]] = []
    chosen_ids = set()
    for sc, cand, evs in designed:
        if len(selected) >= sd.snp_top_n:
            break
        chosen = next((ev for ev in evs if ev.passed), None)
        if chosen is None:
            if any(r.startswith("background") for ev in evs for r in ev.reasons):
                locus_records[sc.locus]["Outcome"] = "all passing pair(s) rejected by the background screen"
            continue
        # Overlapping windows often hold the same amplicon; compare amplicons on the reference.
        amp = _reference_amplicon(chosen, reference.genome_id)
        overlapping = [s.locus for s, _, e in selected if _regions_overlap(amp, _reference_amplicon(e, reference.genome_id))]
        if overlapping:
            locus_records[sc.locus]["Outcome"] = f"skipped: amplicon overlaps selected locus {overlapping[0]}"
            continue
        locus_records[sc.locus]["Outcome"] = f"selected pair {chosen.pair.pair_index}"
        selected.append((sc, cand, chosen))
        chosen_ids.add(id(chosen))
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

    screen_summary = amplicon_screen(selected, genomes, genome_seqs, background_files, primer_dir, cfg) if selected else None

    pair_records = [_pair_record(ev, sc, id(ev) in chosen_ids) for sc, _, evs in designed for ev in evs]
    pd.DataFrame(pair_records).to_csv(primer_dir / "primer_design_summary.csv", index=False)
    pd.DataFrame(list(locus_records.values())).to_csv(primer_dir / "locus_design_log.csv", index=False)
    pd.DataFrame([r for r in pair_records if r["Selected"]]).to_csv(primer_dir / "selected_primers.csv", index=False)

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
        "N_Loci_Designed": len(designed),
        "N_Selected_Primer_Sets": len(selected),
        "Seconds": round(elapsed),
    }


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
    df = pd.DataFrame(rows)
    df.to_csv(cfg.output_dir / "run_summary.csv", index=False)
    return df
