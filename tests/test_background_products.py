"""Distinguishing background products from the target amplicons by sequence."""

import random

from pristine.pipeline import MANY, min_differences
from pristine.seqio import revcomp


def _amplicon(seed: int = 1, n: int = 1200) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def _substitute(seq: str, positions) -> str:
    s = list(seq)
    for i in positions:
        s[i] = {"A": "C", "C": "G", "G": "T", "T": "A"}[s[i]]
    return "".join(s)


def test_identical_and_reverse_complement_are_zero():
    a = _amplicon()
    assert min_differences(a, [a], cap=3) == 0
    assert min_differences(revcomp(a), [a], cap=3) == 0


def test_substitutions_are_counted():
    a = _amplicon()
    assert min_differences(_substitute(a, [100, 500]), [a], cap=10) == 2
    assert min_differences(_substitute(a, [100, 300, 500, 700, 900]), [a], cap=10) == 5


def test_an_indel_run_is_one_event():
    a = _amplicon()
    assert min_differences(a[:600] + a[610:], [a], cap=10) == 1
    assert min_differences(a[:600] + "ACGTACGTAC" + a[600:], [a], cap=10) == 1


def test_closest_target_amplicon_counts():
    a = _amplicon()
    b = _substitute(a, [50, 150, 250, 350])
    assert min_differences(_substitute(a, [50]), [b, a], cap=10) == 1


def test_unrelated_sequence_is_not_aligned():
    assert min_differences(_amplicon(2), [_amplicon(1)], cap=3) == MANY


def test_amplicon_alignment_plot_and_difference_table(tmp_path):
    import pandas as pd
    from pristine.plots import plot_amplicon_alignment

    ref = _amplicon(n=600)
    rows = [("T_1", "target", ref), ("T_2", "target", ref),
            ("N_1", "non-target", _substitute(ref, [10, 20, 30])),
            ("B_1", "background", revcomp(ref[:300] + ref[305:]))]  # a 5 bp deletion, other strand
    png, tsv = tmp_path / "a.png", tmp_path / "a.tsv"
    plot_amplicon_alignment("L", ("T_1", ref), rows, png, tsv)
    assert png.stat().st_size > 0
    diffs = pd.read_csv(tsv, sep="\t")
    assert sorted(diffs.loc[diffs.Sequences == "N_1", "Reference_Position"]) == [11, 21, 31]
    deleted = diffs[(diffs.Sequences == "B_1") & (diffs.Base == "-")]
    assert len(deleted) == 5 and diffs[diffs.Sequences == "T_1;T_2"].empty
