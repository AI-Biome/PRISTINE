import logging
from pathlib import Path

import pytest

from pristine.config import load_config, parse_product_size_range


def test_parse_product_size_range_formats():
    assert parse_product_size_range("500-1800") == [[500, 1800]]
    assert parse_product_size_range("500-1800 2000-2500") == [[500, 1800], [2000, 2500]]
    assert parse_product_size_range([500, 1800]) == [[500, 1800]]
    assert parse_product_size_range([[500, 1800]]) == [[500, 1800]]
    with pytest.raises(ValueError):
        parse_product_size_range("1800-500")


def test_relative_paths_resolve_against_config_dir(tmp_path, monkeypatch):
    (tmp_path / "data").mkdir()
    (tmp_path / "config.yaml").write_text("input_dir: data\noutput_dir: out\nmax_cores: 1\n")
    monkeypatch.chdir(tmp_path.parent)  # cwd deliberately different
    cfg = load_config(tmp_path / "config.yaml")
    assert cfg.input_dir == (tmp_path / "data").resolve()
    assert cfg.output_dir == (tmp_path / "out").resolve()
    assert cfg.primer3_params["PRIMER_MAX_NS_ACCEPTED"] == 0
    assert cfg.product_size_ranges() == [[500, 1800]]
    assert cfg.max_product_size == 1800


def test_legacy_config_keys_are_accepted_with_warning(tmp_path, caplog):
    (tmp_path / "raw").mkdir()
    (tmp_path / "config.yaml").write_text(
        """
input_type: raw
input_paths:
  raw_dir: raw
output_dir: out
max_cores: 8
aligner: mafft
snp_avg_prop_threshold: 0.05
primer3:
  global_params:
    PRIMER_MAX_NS_ACCEPTED: 60
    PRIMER_MIN_TM: 43.0
    PRIMER_MAX_TM: 50.0
    PRIMER_PRODUCT_SIZE_RANGE: "500-1800"
validation:
  perform: no
"""
    )
    with caplog.at_level(logging.WARNING):
        cfg = load_config(tmp_path / "config.yaml")
    assert cfg.input_dir == (tmp_path / "raw").resolve()
    assert cfg.primer3_params["PRIMER_MAX_NS_ACCEPTED"] == 0
    # OPT_TM default (60) lies outside the user's 43-50 range and must be pulled inside
    assert 43.0 <= cfg.primer3_params["PRIMER_OPT_TM"] <= 50.0
    assert any("no longer used" in r.message for r in caplog.records)


def test_window_must_cover_product_size(tmp_path):
    (tmp_path / "raw").mkdir()
    (tmp_path / "config.yaml").write_text(
        "input_dir: raw\noutput_dir: out\nlocus_discovery:\n  window_size: 1000\n  step: 500\n"
    )
    with pytest.raises(ValueError, match="window_size"):
        load_config(tmp_path / "config.yaml")


def test_unknown_key_rejected(tmp_path):
    (tmp_path / "raw").mkdir()
    (tmp_path / "config.yaml").write_text("input_dir: raw\noutput_dir: out\nsnp_primer_design:\n  bogus: 1\n")
    with pytest.raises(ValueError, match="bogus"):
        load_config(tmp_path / "config.yaml")


def test_binding_thresholds_derived_from_primer3_tm(tmp_path):
    (tmp_path / "raw").mkdir()
    (tmp_path / "config.yaml").write_text("input_dir: raw\noutput_dir: out\n")
    cfg = load_config(tmp_path / "config.yaml")
    assert cfg.annealing_temp == 54.0  # PRIMER_MIN_TM 57 - 3
    assert cfg.min_site_tm == 49.0  # minus 5 C margin
    (tmp_path / "config.yaml").write_text(
        "input_dir: raw\noutput_dir: out\nprimer_selection:\n  annealing_temp: 58\n  binding_tm_margin: 3\n"
    )
    assert load_config(tmp_path / "config.yaml").min_site_tm == 55.0
    (tmp_path / "config.yaml").write_text("input_dir: raw\noutput_dir: out\nprimer_selection:\n  use_thermodynamics: false\n")
    assert load_config(tmp_path / "config.yaml").min_site_tm is None


def test_seed_must_cover_three_prime_mismatches(tmp_path):
    """With m mismatches allowed in an R-base 3' region, an exact seed of k bases is
    guaranteed only if R >= (m + 1) * k."""
    (tmp_path / "raw").mkdir()
    (tmp_path / "config.yaml").write_text(
        "input_dir: raw\noutput_dir: out\nprimer_selection:\n  three_prime_region: 8\n  seed_length: 5\n"
    )
    with pytest.raises(ValueError, match="seed_length"):
        load_config(tmp_path / "config.yaml")
    (tmp_path / "config.yaml").write_text(
        "input_dir: raw\noutput_dir: out\nprimer_selection:\n  three_prime_region: 12\n"
        "  three_prime_max_mismatches: 2\n  seed_length: 4\n"
    )
    assert load_config(tmp_path / "config.yaml").primer_selection.seed_length == 4
