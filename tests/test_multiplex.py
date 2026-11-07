import random

from pristine.insilico_pcr import BindingParams, GenomeIndex
from pristine.multiplex import conflicting_loci, cross_pair_products
from pristine.seqio import revcomp

A_L, A_R = "ACGTTGCAGGCTAGCTTACG", "GGATCCATTGACCGTAGCAA"
B_L, B_R = "TTGACCGGATAGCTAGGCAT", "CAGTCGATGCATTGCAGTCA"
STRICT = BindingParams(max_mismatches=0, three_prime_max_mismatches=0)


def _filler(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def _genome(gap_between_loci, seed=3):
    """Locus A (A_L ... A_R) followed by locus B (B_L ... B_R), `gap_between_loci` bp apart."""
    rng = random.Random(seed)
    locus_a = A_L + _filler(rng, 600) + revcomp(A_R)
    locus_b = B_L + _filler(rng, 600) + revcomp(B_R)
    return _filler(rng, 300) + locus_a + _filler(rng, gap_between_loci) + locus_b + _filler(rng, 300)


PAIRS = [("locA", A_L, A_R), ("locB", B_L, B_R)]


def test_close_loci_give_cross_pair_products():
    idx = {"g1": GenomeIndex({"c": _genome(200)})}
    prods = cross_pair_products(PAIRS, idx, STRICT, max_product_size=4000)
    # A_L with B_R spans both loci; B_L with A_R would point away from each other
    assert len(prods) == 1
    row = prods.iloc[0]
    assert (row.Forward_Primer, row.Reverse_Primer) == ("locA_LEFT", "locB_RIGHT")
    assert row.Size == 2 * (20 + 600 + 20) + 200
    assert conflicting_loci(prods) == [("locA", "locB")]


def test_distant_loci_and_single_pairs_give_nothing():
    idx = {"g1": GenomeIndex({"c": _genome(5000)})}
    assert cross_pair_products(PAIRS, idx, STRICT, max_product_size=4000).empty
    # products inside one pair are the per-pair screen's business, not this one's
    assert cross_pair_products(PAIRS[:1], {"g1": GenomeIndex({"c": _genome(200)})}, STRICT, 4000).empty
    assert conflicting_loci(cross_pair_products(PAIRS, idx, STRICT, 4000)) == []
