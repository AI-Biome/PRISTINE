"""Discover input genomes and assign them to target groups.

Input layout (unchanged from the original README)::

    input/
    ├── Genus_species/              # one folder per target group
    │   ├── Genus_species_1.fasta   # target genomes
    │   ├── Genus_species_2.fasta
    │   └── non-targets/            # non-target genomes for this group only
    │       └── Other_species_1.fasta
    └── non-targets/                # non-target genomes used for every group
        └── Third_species_1.fasta

Input files are never modified. Every genome gets a stable ``genome_id`` (its file
stem) and a ``species`` label (the stem without a trailing ``_<number>``).
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

FASTA_EXTS = {".fasta", ".fas", ".fa", ".fna"}
NON_TARGET_DIRNAME = "non-targets"
ROLE_TARGET = "target"
ROLE_GROUP_NONTARGET = "group_nontarget"
ROLE_GLOBAL_NONTARGET = "global_nontarget"


@dataclass(frozen=True)
class Genome:
    genome_id: str
    path: Path
    species: str
    role: str  # one of ROLE_*
    group: str  # target group this genome is used for ("" for global non-targets)

    @property
    def is_target(self) -> bool:
        return self.role == ROLE_TARGET


def species_from_stem(stem: str) -> str:
    """'Yersinia_pestis_12' -> 'Yersinia_pestis'; 'Yersinia_pestis' -> unchanged."""
    return re.sub(r"_\d+$", "", stem)


def is_fasta(path: Path) -> bool:
    return path.is_file() and path.suffix.lower() in FASTA_EXTS


def _list_fasta(folder: Path) -> List[Path]:
    if not folder.is_dir():
        return []
    return sorted(p for p in folder.iterdir() if is_fasta(p))


def discover_groups(input_dir: Path) -> Dict[str, List[Genome]]:
    """Return {group_name: [genomes used for that group]}.

    Each list contains the group's target genomes, its own non-targets and the global
    non-targets. Raises if a group has no targets or if two genomes share an id.
    """
    input_dir = Path(input_dir)
    if not input_dir.is_dir():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")

    global_nt = [
        Genome(p.stem, p.resolve(), species_from_stem(p.stem), ROLE_GLOBAL_NONTARGET, "")
        for p in _list_fasta(input_dir / NON_TARGET_DIRNAME)
    ]

    groups: Dict[str, List[Genome]] = {}
    for folder in sorted(p for p in input_dir.iterdir() if p.is_dir()):
        if folder.name == NON_TARGET_DIRNAME:
            continue
        group = folder.name
        targets = [
            Genome(p.stem, p.resolve(), species_from_stem(p.stem), ROLE_TARGET, group)
            for p in _list_fasta(folder)
        ]
        if not targets:
            raise ValueError(f"Group folder {folder} contains no target FASTA files")
        group_nt = [
            Genome(p.stem, p.resolve(), species_from_stem(p.stem), ROLE_GROUP_NONTARGET, group)
            for p in _list_fasta(folder / NON_TARGET_DIRNAME)
        ]
        genomes = targets + group_nt + [
            Genome(g.genome_id, g.path, g.species, g.role, group) for g in global_nt
        ]
        seen: Dict[str, Path] = {}
        for g in genomes:
            if g.genome_id in seen and seen[g.genome_id] != g.path:
                raise ValueError(
                    f"Duplicate genome id '{g.genome_id}' in group {group}: {seen[g.genome_id]} and {g.path}"
                )
            seen[g.genome_id] = g.path
        # the same file given as group and global non-target: keep one copy
        unique = {}
        for g in genomes:
            unique.setdefault(g.genome_id, g)
        groups[group] = list(unique.values())

    if not groups:
        raise ValueError(f"No target group folders found in {input_dir}")
    return groups


def choose_reference(genomes: List[Genome]) -> Genome:
    """Reference for window tiling: the largest target genome file (most sequence)."""
    targets = [g for g in genomes if g.is_target]
    return max(targets, key=lambda g: g.path.stat().st_size)
