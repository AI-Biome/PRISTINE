"""Optional plots: informativeness heatmap and SNP position histograms."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

from .scoring import LocusScore  # noqa: E402

log = logging.getLogger(__name__)


def plot_informativeness_heatmap(summary_csv: Path, out_png: Path, max_loci: int = 100) -> Path | None:
    df = pd.read_csv(summary_csv)
    prop_cols = [c for c in df.columns if c.startswith("Prop_External_")]
    if df.empty or not prop_cols:
        log.info("No external proportion columns; skipping heatmap")
        return None
    df = df.sort_values("Combined_Informativeness_Score", ascending=False).head(max_loci)
    data = df.set_index("Locus")[prop_cols]
    data.columns = [c.replace("Prop_External_", "") for c in prop_cols]
    fig_h = min(60, max(4, 0.25 * len(data)))
    fig, ax = plt.subplots(figsize=(max(6, 1.2 * len(prop_cols) + 4), fig_h))
    im = ax.imshow(data.values, aspect="auto", cmap="YlOrRd")
    ax.set_xticks(range(len(data.columns)))
    ax.set_xticklabels(data.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(data)))
    ax.set_yticklabels(data.index, fontsize=6)
    ax.set_xlabel("Non-target species")
    ax.set_ylabel("Locus")
    ax.set_title("Proportion of informative sites per locus and non-target species")
    fig.colorbar(im, ax=ax, label="Proportion of informative sites")
    fig.tight_layout()
    fig.savefig(out_png, dpi=120)
    plt.close(fig)
    return out_png


def plot_snp_positions(scores: Sequence[LocusScore], out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for sc in scores:
        fig, ax = plt.subplots(figsize=(12, 3))
        for key, positions in sc.positions_by_key().items():
            if positions:
                ax.hist(positions, bins=min(50, max(1, sc.aln_len // 40)), histtype="step", label=key)
        ax.set_title(f"Informative site positions - {sc.locus}")
        ax.set_xlabel("Alignment position")
        ax.set_ylabel("Count")
        ax.set_xlim(0, sc.aln_len)
        if sc.positions_by_key():
            ax.legend(fontsize="x-small", loc="upper right")
        fig.tight_layout()
        fig.savefig(out_dir / f"snp_positions_{sc.locus}.png", dpi=100)
        plt.close(fig)
