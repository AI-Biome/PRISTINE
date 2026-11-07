"""Multiple sequence alignment of candidate loci with mafft, bounded parallelism."""

from __future__ import annotations

import logging
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

log = logging.getLogger(__name__)


def mafft_available() -> bool:
    return shutil.which("mafft") is not None


def align_one(in_fa: Path, out_fa: Path) -> Optional[Path]:
    cmd = ["mafft", "--auto", "--quiet", "--thread", "1", str(in_fa)]
    try:
        with open(out_fa, "w") as out:
            subprocess.run(cmd, stdout=out, stderr=subprocess.DEVNULL, check=True)
    except subprocess.CalledProcessError:
        log.warning("mafft failed for %s", in_fa.name)
        out_fa.unlink(missing_ok=True)
        return None
    return out_fa


def align_loci(inputs: Sequence[Path], out_dir: Path, threads: int) -> List[Tuple[Path, Optional[Path]]]:
    """Align each FASTA in `inputs`; mafft runs single-threaded, `threads` loci at a time."""
    if not mafft_available():
        raise RuntimeError("mafft not found in PATH (run inside the pixi environment)")
    out_dir.mkdir(parents=True, exist_ok=True)
    jobs = [(p, out_dir / p.name) for p in inputs]
    with ThreadPoolExecutor(max_workers=max(1, threads)) as pool:
        results = list(pool.map(lambda j: align_one(*j), jobs))
    return [(inp, res) for (inp, _), res in zip(jobs, results)]
