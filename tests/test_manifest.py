import hashlib
from pathlib import Path

import pytest

from pristine.manifest import choose_reference, discover_groups, species_from_stem


def test_species_from_stem():
    assert species_from_stem("Yersinia_pestis_12") == "Yersinia_pestis"
    assert species_from_stem("Yersinia_pestis") == "Yersinia_pestis"
    assert species_from_stem("Escherichia_coli_O157_3") == "Escherichia_coli_O157"


def test_discover_groups(dataset):
    root, files = dataset
    groups = discover_groups(root)
    assert list(groups) == ["Targetus_primus"]
    genomes = groups["Targetus_primus"]
    assert sorted(g.genome_id for g in genomes if g.is_target) == ["Targetus_primus_1", "Targetus_primus_2"]
    roles = {g.genome_id: g.role for g in genomes}
    assert roles["Targetus_secundus_1"] == "group_nontarget"
    assert roles["Alienus_tertius_1"] == "global_nontarget"
    assert {g.species for g in genomes} == {"Targetus_primus", "Targetus_secundus", "Alienus_tertius"}
    assert choose_reference(genomes).is_target


def test_discover_groups_does_not_touch_inputs(dataset):
    root, files = dataset
    before = {p: hashlib.md5(p.read_bytes()).hexdigest() for ps in files.values() for p in ps}
    discover_groups(root)
    after = {p: hashlib.md5(p.read_bytes()).hexdigest() for ps in files.values() for p in ps}
    assert before == after


def test_group_without_targets_is_an_error(tmp_path):
    (tmp_path / "Empty_group" / "non-targets").mkdir(parents=True)
    (tmp_path / "Empty_group" / "non-targets" / "X_y_1.fasta").write_text(">a\nACGT\n")
    with pytest.raises(ValueError, match="no target"):
        discover_groups(tmp_path)


def test_no_groups_is_an_error(tmp_path):
    with pytest.raises(ValueError):
        discover_groups(tmp_path)
