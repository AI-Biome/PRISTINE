from pathlib import Path

from pristine.config import LocusDiscoveryParams
from pristine.loci import PafHit, classify_hits, locus_id, parse_paf, sanitize_contig, select_candidate_loci, tile_genome, Window
from pristine.manifest import Genome


def test_tile_genome_covers_tail_and_skips_n_rich():
    g = {"c1": "ACGT" * 600, "c2": "N" * 1500, "short": "ACGT" * 10}
    w = tile_genome(g, window_size=1000, step=400, max_n_fraction=0.05)
    ids = [x.window_id for x in w]
    assert ids[0] == "c1:0-1000"
    assert ids[-1] == "c1:1400-2400"  # tail window appended
    assert not any(x.contig in ("c2", "short") for x in w)


def _hit(qcov=1.0, primary=True, de=0.01, tname="chr", tstart=100):
    qlen = 1000
    return PafHit("w", qlen, 0, int(qlen * qcov), "+", tname, 50000, tstart, tstart + int(qlen * qcov),
                  900, 1000, 60, primary, de)


def test_classify_hits_rules():
    p = LocusDiscoveryParams()
    assert classify_hits([_hit()], p) is not None
    assert classify_hits([], p) is None
    assert classify_hits([_hit(qcov=0.5)], p) is None  # partial
    assert classify_hits([_hit(), _hit(primary=False, tstart=20000)], p) is None  # two copies
    assert classify_hits([_hit(), _hit(qcov=0.2, primary=False, tstart=20000)], p) is not None  # tiny secondary ok
    assert classify_hits([_hit(de=0.3)], p) is None  # too divergent


def test_parse_paf(tmp_path):
    paf = tmp_path / "x.paf"
    paf.write_text("w1\t1000\t0\t1000\t-\tchr\t50000\t100\t1100\t950\t1000\t60\ttp:A:P\tde:f:0.02\n"
                   "w1\t1000\t0\t400\t+\tchr\t50000\t9000\t9400\t380\t400\t0\ttp:A:S\n")
    hits = parse_paf(paf)
    assert len(hits["w1"]) == 2
    assert hits["w1"][0].strand == "-" and hits["w1"][0].divergence == 0.02 and hits["w1"][0].primary
    assert not hits["w1"][1].primary and abs(hits["w1"][1].divergence - 0.05) < 1e-9


def test_select_candidate_loci_requires_presence_and_ranks_by_divergence():
    p = LocusDiscoveryParams(max_candidate_loci=10)
    genomes = [Genome("t1", Path("x"), "T", "target", "G"), Genome("n1", Path("y"), "N", "group_nontarget", "G")]
    windows = [Window("w1", "c", 0, 1000), Window("w2", "c", 500, 1500), Window("w3", "c", 1000, 2000)]
    paf = {
        "t1": {"w1": [_hit(de=0.0)], "w2": [_hit(de=0.0)], "w3": [_hit(de=0.0)]},
        "n1": {"w1": [_hit(de=0.01)], "w2": [_hit(de=0.05)]},  # w3 absent in the non-target
    }
    kept = select_candidate_loci(windows, paf, genomes, p)
    assert [c.window.window_id for c in kept] == ["w2", "w1"]
    p2 = LocusDiscoveryParams(require_in_nontargets=False)
    kept = select_candidate_loci(windows, paf, genomes, p2)
    assert {c.window.window_id for c in kept} == {"w1", "w2", "w3"}


def test_locus_id_is_filesystem_safe_and_short():
    w = Window("x", "accn|CP130901 " + "Yersinia_" * 30, 1000, 3000)
    lid = locus_id(w)
    assert lid.startswith("accn_CP130901_Yersinia")
    assert lid.endswith("_1000_3000") and len(lid) < 60
    assert "|" not in lid and " " not in lid
    assert sanitize_contig("") == "contig"
