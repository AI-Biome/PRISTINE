import random

import pytest

from pristine.config import DEFAULT_PRIMER3_PARAMS
from pristine.primers import best_window, choose_target_region, design_primers


def test_best_window():
    assert best_window([1, 2, 3, 100, 101], 10) == (1, 3)
    assert best_window([], 10) == (None, 0)
    assert best_window([5], 10) == (5, 1)


def test_choose_target_region_merges_windows():
    pos = {"External_a": [100, 110, 120, 130], "External_b": [150, 160, 170]}
    start, length = choose_target_region(pos, window_size=100, min_snps=3, max_target_span=1000, template_len=2000)
    assert start == 100 and start + length == 250


def test_choose_target_region_requires_min_snps():
    assert choose_target_region({"a": [1, 2]}, 100, 3, 1000, 2000) is None


def test_choose_target_region_falls_back_when_span_too_wide():
    pos = {"a": [100, 110, 120], "b": [1500, 1510, 1520, 1530]}
    start, length = choose_target_region(pos, window_size=100, min_snps=3, max_target_span=400, template_len=2000)
    assert length <= 400
    assert start == 1500  # densest window over all positions


def test_design_primers_avoids_masked_positions_and_covers_target():
    rng = random.Random(3)
    seq = "".join(rng.choice("ACGT") for _ in range(1500))
    # mask a variable block and sprinkle Ns
    template = seq[:600] + "N" * 100 + seq[700:]
    params = dict(DEFAULT_PRIMER3_PARAMS)
    params["PRIMER_PRODUCT_SIZE_RANGE"] = [[300, 1200]]
    pairs, explain = design_primers("L", template, (600, 100), params)
    assert pairs, explain
    for p in pairs:
        assert "N" not in p.left_seq and "N" not in p.right_seq
        assert p.amplicon_start <= 600 and p.amplicon_end >= 699
        assert template[p.left_start : p.left_start + p.left_len] == p.left_seq
