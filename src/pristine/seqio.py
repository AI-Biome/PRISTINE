"""Small FASTA helpers (no input mutation, sequences kept upper-case)."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, Iterator, Tuple

_COMPLEMENT = str.maketrans("ACGTRYKMBDHVNacgtrykmbdhvn-", "TGCAYRMKVHDBNtgcayrmkvhdbn-")


def read_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    """Yield (id, sequence) pairs; id is the first whitespace-separated token."""
    name = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].split()[0] if line[1:].strip() else ""
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def read_fasta_dict(path: Path) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for name, seq in read_fasta(path):
        if name in out:
            raise ValueError(f"Duplicate record id '{name}' in {path}")
        out[name] = seq.upper()
    return out


def write_fasta(path: Path, records: Iterable[Tuple[str, str]], width: int = 80) -> None:
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i : i + width] + "\n")


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]
