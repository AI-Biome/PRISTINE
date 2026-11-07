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


# ---- amplicon alignment of target, non-target and background sequences ---------------
# Colours: the four bases take the first four slots of a validated categorical palette
# (blue, orange, aqua, yellow; CVD-separated); matches recede to a light hairline grey.
_SURFACE = "#fcfcfb"
_MATCH = "#e1e0d9"
_BASE_COLOURS = {"A": "#2a78d6", "C": "#eb6834", "G": "#1baf7a", "T": "#eda100"}
_GAP = "#52514e"
_OTHER = "#898781"
_INK = "#0b0b0b"
_INK_SECONDARY = "#52514e"
_CODES = {"A": 1, "C": 2, "G": 3, "T": 4}  # 0 match, 5 gap, 6 other base


def _project(ref: str, seq: str):
    """Align `seq` (the better orientation) to `ref` and project it onto reference
    coordinates: per reference position a code (0 match, 1-4 substituted base, 5 deleted,
    6 other base), plus insertions as {reference position: inserted bases}."""
    from Bio.Align import PairwiseAligner
    from .seqio import revcomp

    aligner = PairwiseAligner(mode="global", match_score=1, mismatch_score=-2, open_gap_score=-4, extend_gap_score=-1)
    seq = seq.upper()
    best = max((aligner.align(ref, s)[0] for s in (seq, revcomp(seq))), key=lambda a: a.score)
    codes = [5] * len(ref)
    insertions = {}
    (ta, qa) = best.aligned
    q_seq = best.sequences[1]
    prev_q = 0
    for (ts, te), (qs, qe) in zip(ta, qa):
        if qs > prev_q and ts > 0:
            insertions[ts] = q_seq[prev_q:qs]
        for i, j in zip(range(ts, te), range(qs, qe)):
            r, b = ref[i], q_seq[j]
            codes[i] = 0 if r == b else _CODES.get(b, 6)
        prev_q = qe
    return codes, insertions


def plot_amplicon_alignment(locus: str, reference: tuple, rows: Sequence[tuple], out_png: Path, out_tsv: Path,
                            note: str = "") -> Path:
    """Colour every difference between each amplicon and the reference amplicon.

    `reference` is (name, sequence); `rows` are (name, group, sequence) with group one of
    "target", "non-target", "background". Identical sequences within a group are drawn
    once. A TSV with one line per difference is written next to the figure."""
    import numpy as np
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch

    ref_name, ref = reference[0], reference[1].upper()
    order = {"target": 0, "non-target": 1, "background": 2}
    collapsed = {}
    for name, group, seq in rows:
        key = (group, seq.upper())
        collapsed.setdefault(key, []).append(name)
    entries = sorted(collapsed.items(), key=lambda kv: (order.get(kv[0][0], 3), kv[1][0]))

    matrix, labels, groups, ins_points, tsv = [], [], [], [], []
    for (group, seq), names in entries:
        codes, insertions = _project(ref, seq)
        row = len(matrix)
        matrix.append(codes)
        labels.append(names[0] + (f"  (+{len(names) - 1} identical)" if len(names) > 1 else ""))
        groups.append(group)
        for pos, bases in insertions.items():
            ins_points.append((pos - 0.5, row))
            tsv.append((";".join(names), group, pos, "", f"insertion {bases}"))
        for i, c in enumerate(codes):
            if c:
                tsv.append((";".join(names), group, i + 1, ref[i], "-" if c == 5 else ("ACGT"[c - 1] if c < 5 else "other")))
    pd.DataFrame(tsv, columns=["Sequences", "Group", "Reference_Position", "Reference_Base", "Base"]).to_csv(
        out_tsv, sep="\t", index=False)

    data = np.array(matrix)
    n_rows = len(matrix)
    fig_h = 1.6 + 0.22 * n_rows
    fig, ax = plt.subplots(figsize=(14, fig_h), facecolor=_SURFACE)
    ax.set_facecolor(_SURFACE)
    cmap = ListedColormap([_MATCH, *_BASE_COLOURS.values(), _GAP, _OTHER])
    ax.imshow(data, aspect="auto", cmap=cmap, vmin=-0.5, vmax=6.5, interpolation="nearest")
    if ins_points:
        xs, ys = zip(*ins_points)
        ax.scatter(xs, ys, marker="|", s=60, color=_INK, linewidths=1.2, zorder=3)

    # group bands: a surface-coloured gap between groups and the group name on the right
    for g in ("target", "non-target", "background"):
        idx = [i for i, gg in enumerate(groups) if gg == g]
        if not idx:
            continue
        if idx[0] > 0:
            ax.axhline(idx[0] - 0.5, color=_SURFACE, linewidth=3)
        # one line, so the labels of adjacent one-row groups cannot overlap
        n_seqs = sum(len(collapsed[(gg, s)]) for (gg, s), _ in entries if gg == g)
        ax.text(len(ref) + len(ref) * 0.01, (idx[0] + idx[-1]) / 2, f"{g} ({n_seqs})",
                va="center", ha="left", fontsize=8, color=_INK_SECONDARY)

    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(labels, fontsize=7, color=_INK)
    ax.set_xlabel(f"Position in the reference amplicon ({ref_name}, {len(ref)} bp)", color=_INK_SECONDARY)
    ax.tick_params(axis="x", colors=_INK_SECONDARY, labelsize=8)
    ax.tick_params(axis="y", length=0)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color("#c3c2b7")
    ax.set_xlim(-0.5, len(ref) - 0.5)

    handles = [Patch(color=_MATCH, label="same as reference")]
    handles += [Patch(color=c, label=f"{b}") for b, c in _BASE_COLOURS.items()]
    handles += [Patch(color=_GAP, label="deleted"), Patch(color=_OTHER, label="other base"),
                plt.Line2D([], [], marker="|", linestyle="", color=_INK, markersize=9, label="insertion")]
    ax.legend(handles=handles, ncol=len(handles), loc="lower left", bbox_to_anchor=(0, 1.0), frameon=False,
              fontsize=8, labelcolor=_INK, handlelength=1.2, borderaxespad=0.3)
    title = f"{locus}: amplicons aligned to the reference amplicon"
    fig.suptitle(title + (f"\n{note}" if note else ""), x=0.01, ha="left", fontsize=10, color=_INK)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, facecolor=_SURFACE)
    plt.close(fig)
    return out_png
