"""Shared fixtures: synthetic genomes with known structure.

Base genome: random 60 kb sequence. Targets are near-identical copies of the base
(0.1 % random substitutions each). Non-target species carry species-specific
substitutions (2-3 %) so that there are informative sites everywhere, plus a 3 kb
segment that is duplicated in every genome (must be rejected as multi-copy).
"""

from __future__ import annotations

import random
from pathlib import Path
from typing import Dict, List

import os

import pytest

from pristine.seqio import write_fasta

ON_GITHUB_ACTIONS = os.environ.get("GITHUB_ACTIONS") == "true"


def pytest_collection_modifyitems(config, items):
    """The suite needs minimap2, mafft and multi-megabase scratch data from the pixi
    environment; it is not run on GitHub Actions."""
    if not ON_GITHUB_ACTIONS:
        return
    skip = pytest.mark.skip(reason="PRISTINE tests are not run on GitHub Actions; run `pixi run test` locally")
    for item in items:
        item.add_marker(skip)


GENOME_LEN = 60_000
DUP_START, DUP_LEN, DUP_INSERT = 10_000, 3_000, 40_000


def random_seq(rng: random.Random, n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(seq: str, rate: float, rng: random.Random) -> str:
    s = list(seq)
    for i in range(len(s)):
        if rng.random() < rate:
            s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
    return "".join(s)


def make_dataset(root: Path, seed: int = 7) -> Dict[str, List[Path]]:
    rng = random.Random(seed)
    base = random_seq(rng, GENOME_LEN)
    # a repeat present twice in every genome
    base = base[:DUP_INSERT] + base[DUP_START : DUP_START + DUP_LEN] + base[DUP_INSERT:]
    species_seqs = {
        "Targetus_primus": base,
        "Targetus_secundus": mutate(base, 0.02, rng),
        "Alienus_tertius": mutate(base, 0.03, rng),
    }
    group = root / "Targetus_primus"
    (group / "non-targets").mkdir(parents=True)
    (root / "non-targets").mkdir()
    files: Dict[str, List[Path]] = {"targets": [], "group_nt": [], "global_nt": []}
    for i in (1, 2):
        p = group / f"Targetus_primus_{i}.fasta"
        write_fasta(p, [("chr", mutate(species_seqs["Targetus_primus"], 0.001, rng))])
        files["targets"].append(p)
    for i in (1, 2):
        p = group / "non-targets" / f"Targetus_secundus_{i}.fasta"
        write_fasta(p, [("chr", mutate(species_seqs["Targetus_secundus"], 0.001, rng))])
        files["group_nt"].append(p)
    p = root / "non-targets" / "Alienus_tertius_1.fasta"
    write_fasta(p, [("chr", species_seqs["Alienus_tertius"])])
    files["global_nt"].append(p)
    return files


@pytest.fixture
def dataset(tmp_path: Path):
    root = tmp_path / "input"
    root.mkdir()
    files = make_dataset(root)
    return root, files


CONFIG_TEMPLATE = """
input_dir: {input_dir}
output_dir: {output_dir}
max_cores: 2
snp_avg_prop_threshold: 0.005
primer3:
  global_params:
    PRIMER_NUM_RETURN: 5
    PRIMER_PRODUCT_SIZE_RANGE: "500-1200"
locus_discovery:
  window_size: 1400
  step: 700
  max_candidate_loci: 40
  refine_with_mafft: {refine}
snp_primer_design:
  snp_window_size: 250
  snp_top_n: 2
  min_snps: 3
  max_loci_to_try: 20
primer_selection:
  conservation: {conservation}
  background_dir: {background}
make_plots: false
"""


def write_config(path: Path, input_dir: str, output_dir: str = "output", refine: bool = False,
                 conservation: str = "all", background: str | None = None) -> Path:
    path.write_text(CONFIG_TEMPLATE.format(input_dir=input_dir, output_dir=output_dir,
                                           refine=str(refine).lower(), conservation=conservation,
                                           background=background if background else "null"))
    return path


@pytest.fixture
def config_file(tmp_path: Path, dataset):
    root, _ = dataset
    return write_config(tmp_path / "config.yaml", root.name)
