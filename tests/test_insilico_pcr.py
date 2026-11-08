import random

import pytest

from pristine.insilico_pcr import BindingParams, GenomeIndex, Thermo, amplify, classify_products, find_sites
from pristine.seqio import revcomp

FWD = "ACGTTGCAGGCTAGCTTACG"  # 20-mer
REV = "GGATCCATTGACCGTAGCAA"
STRICT = BindingParams(max_mismatches=0, three_prime_max_mismatches=0)
LOOSE = BindingParams(max_mismatches=1, three_prime_max_mismatches=0)


def _genome(core_len=800, rng_seed=1):
    rng = random.Random(rng_seed)
    filler = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    return filler(500) + FWD + filler(core_len) + revcomp(REV) + filler(500)


def _flip(base):
    return {"A": "C", "C": "A", "G": "T", "T": "G"}[base]


def test_find_sites_both_strands_and_mismatch_rules():
    g = _genome()
    plus = find_sites("c", g, FWD, "F", STRICT)
    assert [s.strand for s in plus] == ["+"] and plus[0].start == 500 and plus[0].mismatches == 0
    minus = find_sites("c", g, REV, "R", STRICT)
    assert [s.strand for s in minus] == ["-"] and minus[0].end == 500 + 20 + 800 + 20
    # one mismatch at the 5' end is tolerated with max_mismatches=1 but not 0
    mm_fwd = _flip(FWD[0]) + FWD[1:]
    assert find_sites("c", g, mm_fwd, "F", STRICT) == []
    hits = find_sites("c", g, mm_fwd, "F", LOOSE)
    assert len(hits) == 1 and hits[0].mismatches == 1
    # a mismatch at the 3' terminus is rejected when no 3' mismatch is allowed, found when one is
    mm_3 = FWD[:-1] + _flip(FWD[-1])
    assert find_sites("c", g, mm_3, "F", BindingParams(max_mismatches=3, three_prime_max_mismatches=0)) == []
    hits = find_sites("c", g, mm_3, "F", BindingParams(max_mismatches=3, three_prime_max_mismatches=1))
    assert len(hits) == 1 and hits[0].mismatches == 1
    # two mismatches in the 3' region exceed the allowance even if the total is fine
    mm_33 = FWD[:-3] + _flip(FWD[-3]) + FWD[-2] + _flip(FWD[-1])
    assert find_sites("c", g, mm_33, "F", BindingParams(max_mismatches=6, three_prime_max_mismatches=1)) == []
    # GenomeIndex input gives identical results to a plain string
    assert find_sites("c", GenomeIndex({"c": g}), FWD, "F", STRICT) == plus


def test_lowercase_and_N_in_genome():
    g = _genome().lower()
    assert len(find_sites("c", g, FWD, "F", STRICT)) == 1
    g2 = _genome()[:505] + "N" + _genome()[506:]
    assert find_sites("c", g2, FWD, "F", STRICT) == []
    assert len(find_sites("c", g2, FWD, "F", LOOSE)) == 1


def test_amplify_single_product_and_offtarget_detection():
    g = _genome()
    products = amplify({"c": g}, FWD, REV, max_product_size=3000, params=STRICT)
    assert len(products) == 1
    assert products[0].size == 840 and products[0].fwd_primer == "F" and products[0].rev_primer == "R"
    # a second forward site downstream on the minus strand gives an F/F off-target product
    g2 = g + revcomp(FWD) + "ACGT" * 50
    products = amplify({"c": g2}, FWD, REV, max_product_size=3000, params=STRICT)
    assert len(products) == 2
    res = classify_products("g", products, ("c", 400, 1500))
    assert len(res.on_target) == 1 and len(res.off_target) == 1
    assert {p.rev_primer for p in res.off_target} == {"F"}


def test_amplify_respects_max_product_size():
    g = _genome(core_len=5000)
    assert amplify({"c": g}, FWD, REV, max_product_size=3000, params=STRICT) == []
    assert len(amplify({"c": g}, FWD, REV, max_product_size=6000, params=STRICT)) == 1


def test_swapped_primer_labels_still_amplify():
    g = _genome()
    prods = amplify({"c": g}, REV, FWD, max_product_size=3000, params=STRICT)
    assert len(prods) == 1 and prods[0].size == 840
    assert (prods[0].fwd_primer, prods[0].rev_primer) == ("R", "F")


def test_unrelated_primers_give_nothing():
    g = _genome()
    assert amplify({"c": g}, "TTTTTTTTTTTTTTTTTTTT", "CCCCCCCCCCCCCCCCCCCC", max_product_size=3000, params=STRICT) == []


# ---------------------------------------------------------------------------
# thermodynamic binding model
# ---------------------------------------------------------------------------

PRIMER = "GGTATAGAGCCTGACCGGGA"  # Tm ~60 C under PCR salt conditions
THERMO = BindingParams(min_tm=51.0)  # default cutoff: PRIMER_MIN_TM 59 - 3 - 5, Tm at 250 nM primer


def _embed(site, rng_seed=5):
    rng = random.Random(rng_seed)
    filler = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    return filler(300) + site + filler(300)


def test_thermo_perfect_match_tm_matches_primer_tm():
    t = Thermo()
    assert abs(t.duplex_tm(PRIMER, revcomp(PRIMER)) - t.primer_tm(PRIMER)) < 1.0


def test_thermo_single_5prime_mismatch_binds_but_many_do_not():
    one_mm = _flip(PRIMER[0]) + PRIMER[1:]
    hits = find_sites("c", _embed(one_mm), PRIMER, "F", THERMO)
    assert len(hits) == 1 and hits[0].tm is not None and hits[0].tm > 55
    five_mm = "".join(_flip(b) if i in (0, 2, 4, 6, 8) else b for i, b in enumerate(PRIMER))
    assert find_sites("c", _embed(five_mm), PRIMER, "F", THERMO) == []
    # the same five-mismatch site passes a purely count-based model with max_mismatches=6
    assert len(find_sites("c", _embed(five_mm), PRIMER, "F", BindingParams(max_mismatches=6))) == 1


def test_thermo_mismatch_near_3prime_end_counts_as_binding_when_above_threshold():
    # mismatch 5 bases from the 3' end lowers the Tm to ~51 C: still binding at min_tm 49, not at 55
    mm = PRIMER[:15] + _flip(PRIMER[15]) + PRIMER[16:]
    assert len(find_sites("c", _embed(mm), PRIMER, "F", THERMO)) == 1
    strict = BindingParams(min_tm=55.0, thermo=THERMO.thermo)
    assert find_sites("c", _embed(mm), PRIMER, "F", strict) == []


def test_thermo_detects_template_bulge_on_both_strands():
    bulged = PRIMER[:10] + "T" + PRIMER[10:]  # template has one extra base 10 nt from the primer 3' end
    plus = find_sites("c", _embed(bulged), PRIMER, "F", THERMO)
    assert len(plus) == 1 and plus[0].strand == "+" and plus[0].tm > 49
    minus = find_sites("c", _embed(revcomp(bulged)), PRIMER, "F", THERMO)
    assert len(minus) == 1 and minus[0].strand == "-" and minus[0].tm > 49
    # a count-only model without shifts misses the bulge
    assert find_sites("c", _embed(bulged), PRIMER, "F", BindingParams(max_mismatches=2)) == []


def test_thermo_random_seed_hits_are_rejected():
    rng = random.Random(11)
    g = "".join(rng.choice("ACGT") for _ in range(200_000))
    hits = find_sites("c", g, PRIMER, "F", THERMO)
    assert hits == []
    # ...but the real site is found among them
    g2 = g[:100_000] + PRIMER + g[100_000:]
    hits = find_sites("c", g2, PRIMER, "F", THERMO)
    assert [(h.start, h.strand) for h in hits] == [(100_000, "+")]


def test_kmer_codes_match_naive_scan():
    from pristine.insilico_pcr import kmer_code, kmer_codes, _as_bytes
    s = "ACGTNACGTACGGA"
    codes = kmer_codes(_as_bytes(s), 4)
    assert codes.size == len(s) - 3
    # k-mers touching the N are invalid
    assert list(codes[1:5]) == [-1, -1, -1, -1]
    assert codes[0] == kmer_code("ACGT") and codes[5] == kmer_code("ACGT") and codes[-1] == kmer_code("CGGA")
    assert kmer_code("ACNT") == -2


def test_thermo_terminal_mismatch_site_is_reported():
    """A single mismatch at the 3'-terminal base keeps the duplex Tm high; such sites can prime
    (inefficiently) and must be reported rather than ignored."""
    term = PRIMER[:-1] + _flip(PRIMER[-1])
    hits = find_sites("c", _embed(term), PRIMER, "F", THERMO)
    assert len(hits) == 1 and hits[0].mismatches == 1 and hits[0].tm > 49
    minus = find_sites("c", _embed(revcomp(term)), PRIMER, "F", THERMO)
    assert len(minus) == 1 and minus[0].strand == "-"


def test_minus_strand_site_coordinates_and_seed_offsets():
    g = _embed(revcomp(PRIMER))
    hits = find_sites("c", g, PRIMER, "F", BindingParams(max_mismatches=0, three_prime_max_mismatches=0))
    assert len(hits) == 1 and hits[0].strand == "-" and hits[0].start == 300 and hits[0].end == 320
    assert g[hits[0].start : hits[0].end] == revcomp(PRIMER)


# ---------------------------------------------------------------------------
# near misses: products just below the binding cutoff
# ---------------------------------------------------------------------------

def test_near_miss_product_reported_only_on_request():
    from pristine.insilico_pcr import is_near_miss

    rng = random.Random(9)
    filler = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    five_mm = "".join(_flip(b) if i in (0, 2, 4, 6, 8) else b for i, b in enumerate(PRIMER))
    genome = {"c": filler(300) + five_mm + filler(700) + revcomp(REV) + filler(300)}
    site_tm = Thermo().duplex_tm(PRIMER, revcomp(five_mm))
    assert 46.0 <= site_tm < 51.0  # inside the near-miss band of the default cutoff
    params = BindingParams(min_tm=51.0, near_miss_tm=46.0)
    assert amplify(genome, PRIMER, REV, 2000, params) == []
    near = amplify(genome, PRIMER, REV, 2000, params, include_near_misses=True)
    assert len(near) == 1 and is_near_miss(near[0], params)
    assert near[0].weakest_tm == pytest.approx(min(near[0].fwd_tm, near[0].rev_tm))
    # a perfect product is returned either way and is not a near miss
    perfect = {"c": filler(300) + PRIMER + filler(700) + revcomp(REV) + filler(300)}
    prods = amplify(perfect, PRIMER, REV, 2000, params, include_near_misses=True)
    assert len(prods) == 1 and not is_near_miss(prods[0], params)
    assert amplify(perfect, PRIMER, REV, 2000, params) == prods
    # below the near-miss band nothing is reported
    assert amplify(genome, PRIMER, REV, 2000, BindingParams(min_tm=56.0, near_miss_tm=54.0),
                   include_near_misses=True) == []


def test_near_miss_cutoff_must_not_exceed_binding_cutoff():
    with pytest.raises(ValueError, match="near_miss_tm"):
        BindingParams(min_tm=50.0, near_miss_tm=52.0)
    with pytest.raises(ValueError, match="near_miss_tm"):
        BindingParams(near_miss_tm=45.0)


def test_default_primer_concentration_raises_tm():
    # 250 nM (the PRISTINE default) vs primer3's 50 nM: about 2 C higher
    assert 1.5 < Thermo().primer_tm(PRIMER) - Thermo(dna_conc=50.0).primer_tm(PRIMER) < 3.0
