"""End-to-end run on synthetic genomes with the real minimap2 / mafft / primer3."""

import hashlib
import shutil

import pandas as pd
import pytest

from pristine.config import load_config
from pristine.insilico_pcr import BindingParams, amplify
from pristine.pipeline import run_pipeline
from pristine.seqio import read_fasta_dict

from conftest import write_config

needs_tools = pytest.mark.skipif(
    shutil.which("minimap2") is None or shutil.which("mafft") is None,
    reason="minimap2 and mafft required",
)


@needs_tools
def test_end_to_end(config_file, dataset):
    root, files = dataset
    checksums = {p: hashlib.md5(p.read_bytes()).hexdigest() for ps in files.values() for p in ps}

    cfg = load_config(config_file)
    summary = run_pipeline(cfg)
    assert "Error" not in summary.columns, summary.to_string()
    row = summary.iloc[0]
    assert row["N_Selected_Primer_Sets"] >= 1

    out = cfg.output_dir / "Targetus_primus"
    snp = pd.read_csv(out / "informative_loci" / "snp_summary.csv")
    assert {"Prop_External_Targetus_secundus", "Prop_External_Alienus_tertius"} <= set(snp.columns)
    assert (snp["Abs_External_Alienus_tertius"] > 0).all()
    # the duplicated 3 kb segment must not yield a locus
    for locus in snp["Locus"]:
        _, start, end = locus.rsplit("_", 2)
        assert not (int(start) >= 9_900 and int(end) <= 13_100), f"repeat region kept: {locus}"
        assert not (int(start) >= 39_900 and int(end) <= 43_100), f"repeat copy kept: {locus}"

    sel = pd.read_csv(out / "primers" / "selected_primers.csv")
    assert len(sel) == row["N_Selected_Primer_Sets"]
    # selected loci must not overlap each other on the reference (windows overlap by step)
    spans = sorted((int(l.rsplit("_", 2)[1]), int(l.rsplit("_", 2)[2])) for l in sel["Locus"])
    for (s1, e1), (s2, e2) in zip(spans, spans[1:]):
        assert e1 <= s2, spans
    all_pairs = pd.read_csv(out / "primers" / "primer_design_summary.csv")
    assert all_pairs["Selected"].sum() == len(sel)
    # ranked best first: Rank 1..n, background-free before fallback, then locus score, then penalty
    assert list(sel["Rank"]) == list(range(1, len(sel) + 1))
    key = list(zip(sel["Background_Products"].fillna("").astype(str) != "", -sel["Locus_Combined_Score"],
                   sel["Primer3_Penalty"]))
    assert key == sorted(key)
    assert list(all_pairs["Rank"].head(len(sel))) == list(range(1, len(sel) + 1))
    assert all_pairs["Rank"].iloc[len(sel):].isna().all()
    ranked_fa = list(read_fasta_dict(out / "primers" / "selected_primers.fasta").items())
    assert ranked_fa[0] == (f"rank1_{sel['Locus'].iloc[0]}_LEFT", sel["LEFT_PRIMER"].iloc[0])
    assert [h for h, _ in ranked_fa] == [f"rank{r}_{l}_{side}" for r, l in zip(sel["Rank"], sel["Locus"])
                                         for side in ("LEFT", "RIGHT")]
    for _, r in sel.iterrows():
        assert r["Passed_Checks"]
        assert set(r["LEFT_PRIMER"]) <= set("ACGT") and set(r["RIGHT_PRIMER"]) <= set("ACGT")
        assert 57 <= r["LEFT_TM"] <= 63 and 57 <= r["RIGHT_TM"] <= 63
        assert r["Sites_In_Amplicon_External_Alienus_tertius"] >= 3
        assert r["Sites_In_Amplicon_External_Targetus_secundus"] >= 3
        # independent re-check: exactly one product of sane size in every genome
        for ps in files.values():
            for p in ps:
                prods = amplify(read_fasta_dict(p), r["LEFT_PRIMER"], r["RIGHT_PRIMER"], 4000,
                                BindingParams(min_tm=49.0))
                assert len(prods) == 1, (p, prods)
                assert 450 <= prods[0].size <= 1320
        amp = read_fasta_dict(out / "primers" / "amplicons" / f"{r['Locus']}_amplicons.fasta")
        assert len(amp) == 5
        assert (out / "primers" / f"{r['Locus']}_primers.fasta").exists()

    assert (out / "primers" / "locus_design_log.csv").exists()
    assert (out / "informative_loci" / "consensus_sequences").is_dir()

    # inputs untouched
    assert checksums == {p: hashlib.md5(p.read_bytes()).hexdigest() for p in checksums}


@needs_tools
def test_cli_runs(config_file, tmp_path):
    from pristine.cli import main
    rc = main(["--config", str(config_file), "--threads", "1", "--output", str(tmp_path / "cli_out")])
    assert rc == 0
    assert (tmp_path / "cli_out" / "run_summary.csv").exists()
    assert (tmp_path / "cli_out" / "pristine.log").exists()


@needs_tools
def test_end_to_end_with_mafft_refinement_and_background(dataset, tmp_path):
    """refine_with_mafft re-aligns design loci; a background genome that contains the
    reference verbatim must make every pair fail (it would amplify there too)."""
    root, files = dataset
    bg = tmp_path / "background"
    bg.mkdir()
    shutil.copy(files["targets"][0], bg / "Contaminant_1.fasta")
    cfg = load_config(write_config(tmp_path / "config_bg.yaml", root.name, "out_bg", refine=True, background=bg.name))
    summary = run_pipeline(cfg)
    assert "Error" not in summary.columns, summary.to_string()
    assert summary.iloc[0]["N_Selected_Primer_Sets"] == 0
    pairs = pd.read_csv(cfg.output_dir / "Targetus_primus" / "primers" / "primer_design_summary.csv")
    assert len(pairs) > 0
    assert pairs["Failure_Reasons"].str.contains("background Contaminant_1").any()
    passed_group = pairs[~pairs["Failure_Reasons"].str.contains("no product|off-target|products at locus|informative", regex=True)]
    assert passed_group["Failure_Reasons"].str.contains("background Contaminant_1").all()
    log = pd.read_csv(cfg.output_dir / "Targetus_primus" / "primers" / "locus_design_log.csv")
    assert log["Outcome"].str.contains("background").any()
    assert (cfg.output_dir / "Targetus_primus" / "loci" / "mafft").is_dir()


@needs_tools
def test_conservation_targets_mode(dataset, tmp_path):
    root, _ = dataset
    cfg = load_config(write_config(tmp_path / "config_t.yaml", root.name, "out_t", conservation="targets"))
    summary = run_pipeline(cfg)
    assert "Error" not in summary.columns, summary.to_string()
    assert summary.iloc[0]["N_Selected_Primer_Sets"] >= 1


@needs_tools
def test_amplicon_screen_against_background(dataset, tmp_path):
    """A divergent background genome must not reject the primers but shows up in the amplicon screen."""
    import random
    from pristine.seqio import read_fasta_dict, write_fasta

    root, files = dataset
    bg = tmp_path / "background"
    bg.mkdir()
    rng = random.Random(3)
    far = "".join(rng.choice("ACGT") for _ in range(50_000))
    write_fasta(bg / "Unrelated_1.fasta", [("chr", far)])
    cfg = load_config(write_config(tmp_path / "config_as.yaml", root.name, "out_as", background=bg.name))
    summary = run_pipeline(cfg)
    assert "Error" not in summary.columns, summary.to_string()
    assert summary.iloc[0]["N_Selected_Primer_Sets"] >= 1
    screen = pd.read_csv(cfg.output_dir / "Targetus_primus" / "primers" / "amplicon_screen" / "amplicon_screen_summary.csv")
    assert len(screen) == summary.iloc[0]["N_Selected_Primer_Sets"]
    assert (screen["Tolerated_Identity"] < 100).all()
    assert screen["Background_Max_Identity"].isna().all() or (screen["Background_Max_Identity"] < 90).all()
    assert (cfg.output_dir / "Targetus_primus" / "primers" / "amplicon_screen" / "reference_amplicons.fasta").exists()


def _locus_span(locus: str):
    contig, start, end = locus.rsplit("_", 2)
    return contig, int(start), int(end)


@needs_tools
def test_design_continues_until_snp_top_n_after_background_rejections(dataset, tmp_path):
    """If the background screen rejects every locus of the first design round, further
    loci are designed until snp_top_n loci are selected."""
    from pristine.seqio import write_fasta

    root, files = dataset
    # Run A without background and snp_top_n=2 designs loci until two passing regions are
    # found. Run B (snp_top_n=1, with background) designs the same loci in its first round
    # (it asks for 2 * missing passing regions), so all of them are rejected there.
    cfg = load_config(write_config(tmp_path / "config_a.yaml", root.name, "out_a"))
    cfg.snp_primer_design.snp_top_n = 2
    summary = run_pipeline(cfg)
    assert "Error" not in summary.columns, summary.to_string()
    reference = summary.iloc[0]["Reference"]
    pairs = pd.read_csv(cfg.output_dir / "Targetus_primus" / "primers" / "primer_design_summary.csv")
    passing = sorted(set(pairs.loc[pairs["Passed_Checks"], "Locus"]))
    assert passing

    # Background: the reference sequence of every locus that passed in run A.
    ref = read_fasta_dict(next(p for ps in files.values() for p in ps if p.stem == reference))["chr"]
    bg = tmp_path / "background"
    bg.mkdir()
    write_fasta(bg / "Contaminant_1.fasta",
                [(l, ref[max(0, s - 100): e + 100]) for l, (_, s, e) in ((l, _locus_span(l)) for l in passing)])

    cfg_b = load_config(write_config(tmp_path / "config_b.yaml", root.name, "out_b", background=bg.name))
    cfg_b.snp_primer_design.snp_top_n = 1
    summary_b = run_pipeline(cfg_b)
    assert "Error" not in summary_b.columns, summary_b.to_string()
    row = summary_b.iloc[0]
    assert row["N_Selected_Primer_Sets"] == 1
    sel = pd.read_csv(cfg_b.output_dir / "Targetus_primus" / "primers" / "selected_primers.csv")
    assert sel.iloc[0]["Locus"] not in passing
    log = pd.read_csv(cfg_b.output_dir / "Targetus_primus" / "primers" / "locus_design_log.csv")
    assert log["Outcome"].str.contains("rejected by the background screen").any()
    assert not (cfg_b.output_dir / "Targetus_primus" / "primers" / "amplicon_screen" / "tmp").exists()


@needs_tools
def test_nontarget_lacking_locus_with_require_in_nontargets_false(tmp_path):
    """require_in_nontargets: false with conservation: all must allow loci that a
    non-target genome lacks; that genome then simply must give no product."""
    from conftest import make_dataset
    from pristine.seqio import write_fasta

    root = tmp_path / "input"
    root.mkdir()
    files = make_dataset(root)
    # Alienus_tertius keeps only the first 15 kb: most loci are absent from it.
    nt = files["global_nt"][0]
    write_fasta(nt, [("chr", read_fasta_dict(nt)["chr"][:15_000])])
    cfg = load_config(write_config(tmp_path / "config_r.yaml", root.name, "out_r"))
    cfg.locus_discovery.require_in_nontargets = False
    cfg.locus_discovery.max_candidate_loci = 200
    cfg.snp_primer_design.max_loci_to_try = 100
    cfg.snp_primer_design.snp_top_n = 8
    summary = run_pipeline(cfg)
    assert "Error" not in summary.columns, summary.to_string()
    sel = pd.read_csv(cfg.output_dir / "Targetus_primus" / "primers" / "selected_primers.csv")
    assert (sel["Locus"].map(lambda l: _locus_span(l)[1]) >= 16_000).any(), sel["Locus"].tolist()
    for _, r in sel.iterrows():
        prods = amplify(read_fasta_dict(nt), r["LEFT_PRIMER"], r["RIGHT_PRIMER"], 4000, BindingParams(min_tm=49.0))
        if _locus_span(r["Locus"])[1] >= 16_000:
            assert prods == [], (r["Locus"], prods)


@needs_tools
def test_parallel_run_matches_serial_run(dataset, tmp_path):
    """max_cores only changes speed: the parallel in silico PCR (batched loci, pair-level
    worker processes, parallel background screen) must give the same tables and files as a
    one-core run."""
    from pristine.seqio import write_fasta

    root, files = dataset
    # A background genome holding the best passing locus rejects it, so the background
    # screen has work and design needs a second round.
    cfg0 = load_config(write_config(tmp_path / "config_0.yaml", root.name, "out_0"))
    summary0 = run_pipeline(cfg0)
    reference = summary0.iloc[0]["Reference"]
    sel0 = pd.read_csv(cfg0.output_dir / "Targetus_primus" / "primers" / "selected_primers.csv")
    _, s, e = _locus_span(sel0.iloc[0]["Locus"])
    ref = read_fasta_dict(next(p for ps in files.values() for p in ps if p.stem == reference))["chr"]
    bg = tmp_path / "background"
    bg.mkdir()
    write_fasta(bg / "Contaminant_1.fasta", [("locus", ref[max(0, s - 100): e + 100])])
    write_fasta(bg / "Unrelated_1.fasta", [("chr", ref[::-1])])

    outputs = {}
    for cores in (1, 4):
        cfg = load_config(write_config(tmp_path / f"config_{cores}.yaml", root.name, f"out_{cores}", background=bg.name))
        cfg.max_cores = cores
        summary = run_pipeline(cfg)
        assert "Error" not in summary.columns, summary.to_string()
        out = cfg.output_dir / "Targetus_primus"
        outputs[cores] = {
            "summary": summary.drop(columns=["Seconds"]),
            **{name: pd.read_csv(out / "primers" / name)
               for name in ("primer_design_summary.csv", "locus_design_log.csv", "selected_primers.csv",
                            "multiplex_check.csv")},
            "consensus": sorted(p.name for p in (out / "informative_loci" / "consensus_sequences").iterdir()),
        }
    serial, parallel = outputs[1], outputs[4]
    assert serial["primer_design_summary.csv"]["Failure_Reasons"].fillna("").str.contains("background").any()
    for key in serial:
        if isinstance(serial[key], pd.DataFrame):
            pd.testing.assert_frame_equal(serial[key], parallel[key], check_like=False)
        else:
            assert serial[key] == parallel[key], key


@needs_tools
def test_background_products_distinguishable_by_sequence(dataset, tmp_path):
    """A background genome that amplifies like a non-target (the primer sites are conserved)
    is accepted with background_products: distinguish, because its amplicon carries the
    non-target's informative sites, and fails every such pair with background_products: reject."""
    root, files = dataset
    bg = tmp_path / "background"
    bg.mkdir()
    shutil.copy(files["global_nt"][0], bg / "Bystander_1.fasta")  # a copy of Alienus_tertius

    cfg = load_config(write_config(tmp_path / "config_d.yaml", root.name, "out_d", background=bg.name))
    assert cfg.primer_selection.background_products == "distinguish"
    summary = run_pipeline(cfg)
    assert "Error" not in summary.columns, summary.to_string()
    assert summary.iloc[0]["N_Selected_Primer_Sets"] >= 1
    sel = pd.read_csv(cfg.output_dir / "Targetus_primus" / "primers" / "selected_primers.csv")
    assert sel["Background_Products"].fillna("").str.contains("Bystander_1").all()
    # PRISCREEN exports: background products are panel records, the sites table names them
    from pristine.seqio import read_fasta
    primers = cfg.output_dir / "Targetus_primus" / "primers"
    heads = [h for h, _ in read_fasta(primers / "priscreen_panel.fasta")]
    assert any("|species=Bystander|genome=Bystander_1|role=background" in h for h in heads)
    assert any("role=target" in h for h in heads)
    assert (cfg.output_dir / "priscreen_panel.fasta").exists()
    for locus in sel["Locus"]:
        sites = pd.read_csv(primers / "diagnostic_sites" / f"{locus}.tsv", sep="\t")
        assert "Bystander" in sites.columns and len(sites) >= cfg.snp_primer_design.min_snps

    cfg_r = load_config(write_config(tmp_path / "config_r.yaml", root.name, "out_r", background=bg.name))
    cfg_r.primer_selection.background_products = "reject"
    summary_r = run_pipeline(cfg_r)
    assert summary_r.iloc[0]["N_Selected_Primer_Sets"] == 0
    pairs = pd.read_csv(cfg_r.output_dir / "Targetus_primus" / "primers" / "primer_design_summary.csv")
    assert pairs["Failure_Reasons"].fillna("").str.contains("background Bystander_1").any()
    log = pd.read_csv(cfg.output_dir / "Targetus_primus" / "primers" / "locus_design_log.csv")
    assert log["Outcome"].str.contains("fallback: amplifies in the background").any()


@needs_tools
def test_loci_without_background_products_are_preferred(dataset, tmp_path):
    """A locus that amplifies (distinguishably) in the background is only used when there
    are not enough loci without background products."""
    from pristine.seqio import write_fasta

    root, files = dataset
    cfg0 = load_config(write_config(tmp_path / "config_0.yaml", root.name, "out_0"))
    run_pipeline(cfg0)
    sel0 = pd.read_csv(cfg0.output_dir / "Targetus_primus" / "primers" / "selected_primers.csv")
    top = sel0.iloc[0]["Locus"]
    _, s, e = _locus_span(top)
    # the non-target Alienus_tertius has the same coordinates (substitutions only)
    alien = read_fasta_dict(files["global_nt"][0])["chr"]
    bg = tmp_path / "background"
    bg.mkdir()
    write_fasta(bg / "Bystander_1.fasta", [("locus", alien[max(0, s - 100): e + 100])])

    cfg = load_config(write_config(tmp_path / "config_p.yaml", root.name, "out_p", background=bg.name))
    summary = run_pipeline(cfg)
    assert summary.iloc[0]["N_Selected_Primer_Sets"] == cfg.snp_primer_design.snp_top_n
    sel = pd.read_csv(cfg.output_dir / "Targetus_primus" / "primers" / "selected_primers.csv")
    assert top not in set(sel["Locus"])
    assert sel["Background_Products"].isna().all()
    log = pd.read_csv(cfg.output_dir / "Targetus_primus" / "primers" / "locus_design_log.csv").set_index("Locus")
    assert "kept in reserve" in log.loc[top, "Outcome"]


@needs_tools
def test_near_miss_products_reject_the_pair(dataset, tmp_path):
    """A background genome whose copy of a selected primer site binds just below the cutoff
    gives a near-miss product, which rejects the pair like a product would."""
    import itertools

    from pristine.insilico_pcr import is_near_miss
    from pristine.pipeline import binding_params_from_config
    from pristine.seqio import write_fasta

    root, files = dataset
    cfg0 = load_config(write_config(tmp_path / "config_0.yaml", root.name, "out_0"))
    run_pipeline(cfg0)
    top = pd.read_csv(cfg0.output_dir / "Targetus_primus" / "primers" / "selected_primers.csv").iloc[0]
    left, right = top["LEFT_PRIMER"], top["RIGHT_PRIMER"]
    alien = read_fasta_dict(files["global_nt"][0])["chr"]
    pos = alien.find(left)
    assert pos >= 0, "primer sites are conserved in the non-target"
    binding = binding_params_from_config(cfg0)

    # two mismatches in the 5' half of the left primer site that put it into the near-miss band
    def mutate(seq, i, j):
        s = list(seq)
        for k in (pos + i, pos + j):
            s[k] = {"A": "C", "C": "A", "G": "T", "T": "G"}[s[k]]
        return "".join(s)

    lo, hi = pos - 100, pos + 2100
    bystander = None
    for i, j in itertools.combinations(range(len(left) // 2), 2):
        seq = mutate(alien, i, j)[max(0, lo):hi]
        prods = amplify({"locus": seq}, left, right, 4000, binding, include_near_misses=True)
        if prods and all(is_near_miss(p, binding) for p in prods):
            bystander = seq
            break
    assert bystander is not None
    bg = tmp_path / "background"
    bg.mkdir()
    write_fasta(bg / "Bystander_1.fasta", [("locus", bystander)])

    cfg = load_config(write_config(tmp_path / "config_n.yaml", root.name, "out_n", background=bg.name))
    summary = run_pipeline(cfg)
    assert summary.iloc[0]["N_Selected_Primer_Sets"] == cfg.snp_primer_design.snp_top_n
    primers = cfg.output_dir / "Targetus_primus" / "primers"
    pairs = pd.read_csv(primers / "primer_design_summary.csv")
    hit = pairs[(pairs["LEFT_PRIMER"] == left) & (pairs["RIGHT_PRIMER"] == right)]
    assert len(hit) and not hit["Passed_Checks"].any()
    assert hit["Failure_Reasons"].str.contains("background Bystander_1: 1 near-miss product").all()
    assert "Near_Miss_Products" not in pairs.columns
    sel = pd.read_csv(primers / "selected_primers.csv")
    assert (sel["Recommended_Annealing_Temp"] >= cfg.annealing_temp).all()
    assert not ((sel["LEFT_PRIMER"] == left) & (sel["RIGHT_PRIMER"] == right)).any()


@needs_tools
def test_poolable_selection_has_no_cross_pair_products(dataset, tmp_path):
    """With poolable: true (default) no two selected pairs form a product or near miss
    together; with false, close loci may be selected and are only reported."""
    root, _ = dataset
    cfg = load_config(write_config(tmp_path / "config_p.yaml", root.name, "out_p"))
    assert cfg.primer_selection.poolable
    run_pipeline(cfg)
    primers = cfg.output_dir / "Targetus_primus" / "primers"
    assert pd.read_csv(primers / "multiplex_check.csv").empty

    cfg_f = load_config(write_config(tmp_path / "config_f.yaml", root.name, "out_f"))
    cfg_f.primer_selection.poolable = False
    run_pipeline(cfg_f)
    cross = pd.read_csv(cfg_f.output_dir / "Targetus_primus" / "primers" / "multiplex_check.csv")
    log = pd.read_csv(primers / "locus_design_log.csv")
    skipped = log["Outcome"].fillna("").str.contains("cross-pair products with selected locus")
    # every conflict that the unconstrained selection runs into is avoided by skipping a locus
    assert cross.empty or skipped.any()
