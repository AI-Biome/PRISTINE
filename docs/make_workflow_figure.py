#!/usr/bin/env python3
"""Draw the PRISTINE workflow figure as a standalone SVG.

Usage:  python docs/make_workflow_figure.py [output.svg]   (default: PRISTINE-1.svg in the repo root)
For a PDF:  inkscape PRISTINE-1.svg --export-type=pdf
"""

from __future__ import annotations

import sys
from pathlib import Path

FONT = "Helvetica, Arial, 'Liberation Sans', 'DejaVu Sans', sans-serif"
INK, MUTED, ARROW = "#1a1a1a", "#555555", "#3a3a3a"

# stage colours (light panel, dark accent); muted, colour-blind friendly
STAGES = {
    "A": ("Locus discovery", "#eaf1f8", "#2c6aa0"),
    "B": ("Marker scoring", "#ebf5ee", "#2e7d4f"),
    "C": ("Primer design", "#fbf1e4", "#b86e12"),
    "D": ("Specificity validation", "#f3ecf6", "#7a4a96"),
}
DATA_IN = ("#e9eef7", "#33518a")   # input data
DATA_OUT = ("#e8f4ec", "#2e7d4f")  # results

STEPS = [  # (stage, title, detail)
    ("A", "Reference windows", "largest target genome tiled into amplicon-sized windows"),
    ("A", "Genome-wide mapping", "windows aligned to every genome (minimap2)"),
    ("A", "Single-copy filter", "one full-length hit per genome, no repeats"),
    ("B", "Locus alignment", "built from minimap2 cs tags (mafft optional)"),
    ("B", "Informative-site scoring", "substitutions and indel events vs. each non-target"),
    ("C", "Primer design", "conserved sites flanking the SNP-dense region (Primer3)"),
    ("D", "In silico PCR", "all input genomes; 3'-seeded sites, duplex Tm"),
    ("D", "Background screen", "a product too similar to the target amplicons rejects a pair"),
    ("D", "Locus selection", "one marker per region; background-amplifying loci last"),
    ("D", "Multiplex check", "products formed by primers of two selected pairs"),
    ("D", "Amplicon similarity screen", "identity to background genomes and NCBI nt"),
]
PARALLEL = {1, 6, 7}   # steps that use max_cores: concurrent minimap2 jobs, in silico PCR worker processes
I_DESIGN, I_SELECT = 5, 8   # steps 5-7 repeat in rounds until snp_top_n loci are selected

W, H = 720, 1080
CX, BW, BH = 270, 380, 46          # main column centre, step box size
ROW, GAP = 62, 48                  # row pitch within a stage, extra gap between stages
PANEL_X, PANEL_W = 30, 480         # stage panel extent
PANEL_HEAD, PANEL_PAD = 30, 12     # header height and bottom padding of a stage panel
SIDE_X = 610                       # centre of the right-hand column

svg: list[str] = []
add = svg.append


def esc(s: str) -> str:
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def text(x, y, s, size=12, color=INK, weight="normal", anchor="middle", extra=""):
    add(f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" fill="{color}" font-weight="{weight}" '
        f'text-anchor="{anchor}"{extra}>{esc(s)}</text>')


def arrow(points, dashed=False, color=ARROW):
    d = " ".join(("M" if i == 0 else "L") + f"{x:.1f},{y:.1f}" for i, (x, y) in enumerate(points))
    dash = ' stroke-dasharray="4,3"' if dashed else ""
    add(f'<path d="{d}" fill="none" stroke="{color}" stroke-width="1.3"{dash} marker-end="url(#head)"/>')


def document(cx, cy, w, h, title, detail, colors, dashed=False, stack=False):
    """Page symbol with a folded corner; `stack` draws two pages behind it."""
    fill, stroke = colors
    dash = ' stroke-dasharray="4,3"' if dashed else ""
    x, y, f = cx - w / 2, cy - h / 2, 11

    def page(px, py, sw):
        add(f'<path d="M{px},{py} H{px + w - f} L{px + w},{py + f} V{py + h} H{px} Z" fill="{fill}" '
            f'stroke="{stroke}" stroke-width="{sw}"{dash}/>')

    if stack:
        for off in (8, 4):
            page(x + off, y + off, 1.0)
    page(x, y, 1.3)
    add(f'<path d="M{x + w - f},{y} V{y + f} H{x + w}" fill="none" stroke="{stroke}" stroke-width="1"/>')
    text(cx, cy - 2, title, 12.5, INK, "bold")
    text(cx, cy + 13, detail, 10, MUTED)


def database(cx, cy, w, h, title, detail, colors, dashed=False):
    fill, stroke = colors
    dash = ' stroke-dasharray="4,3"' if dashed else ""
    x, y, ry = cx - w / 2, cy - h / 2, 7
    add(f'<path d="M{x},{y + ry} V{y + h - ry} A{w / 2},{ry} 0 0 0 {x + w},{y + h - ry} V{y + ry}" '
        f'fill="{fill}" stroke="{stroke}" stroke-width="1.3"{dash}/>')
    add(f'<ellipse cx="{cx}" cy="{y + ry}" rx="{w / 2}" ry="{ry}" fill="{fill}" stroke="{stroke}" stroke-width="1.3"{dash}/>')
    text(cx, cy + 4, title, 12.5, INK, "bold")
    text(cx, cy + 18, detail, 10, MUTED)


def tag(x_right, cy, color, label="max_cores"):
    """Small pill, right-aligned at x_right: the step runs in parallel worker processes."""
    w = 8 + 6.0 * len(label) + 8
    add(f'<rect x="{x_right - w:.1f}" y="{cy - 7:.1f}" width="{w:.1f}" height="14" rx="7" fill="#ffffff" '
        f'stroke="{color}" stroke-width="1"/>')
    text(x_right - w / 2, cy + 3.5, label, 9, color, "bold")


# ---------------------------------------------------------------------------
add(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" '
    f'font-family="{FONT}">')
add("<title>PRISTINE workflow</title>")
add(f'<defs><marker id="head" viewBox="0 0 10 10" refX="9.5" refY="5" markerWidth="6.5" markerHeight="6.5" '
    f'orient="auto"><path d="M0,0.8 L10,5 L0,9.2 Z" fill="{ARROW}"/></marker></defs>')
add(f'<rect width="{W}" height="{H}" fill="#ffffff"/>')

# inputs
in_y = 52
document(CX - 92, in_y, 160, 46, "Target genomes", "one folder per group", DATA_IN, stack=True)
document(CX + 92, in_y, 160, 46, "Non-target genomes", "close relatives", DATA_IN, stack=True)

# step positions
ys, y, prev = [], 160, None
for stage, _, _ in STEPS:
    if prev is not None:
        y += ROW + (GAP if stage != prev else 0)
    ys.append(y)
    prev = stage

# stage panels (drawn first, behind the steps)
for key, (label, light, dark) in STAGES.items():
    idx = [i for i, s in enumerate(STEPS) if s[0] == key]
    top, bottom = ys[idx[0]] - BH / 2 - PANEL_HEAD, ys[idx[-1]] + BH / 2 + PANEL_PAD
    add(f'<rect x="{PANEL_X}" y="{top:.1f}" width="{PANEL_W}" height="{bottom - top:.1f}" rx="8" fill="{light}"/>')
    add(f'<circle cx="{PANEL_X + 18}" cy="{top + 15:.1f}" r="8.5" fill="none" stroke="{dark}" stroke-width="1.3"/>')
    text(PANEL_X + 18, top + 19, key, 10.5, dark, "bold")
    text(PANEL_X + 33, top + 19, label.upper(), 10.5, dark, "bold", "start", ' letter-spacing="1"')

# arrows from inputs into step 1
arrow([(CX - 70, in_y + 24), (CX - 30, ys[0] - BH / 2 - 2)])
arrow([(CX + 70, in_y + 24), (CX + 30, ys[0] - BH / 2 - 2)])

# steps
for i, (stage, title, detail) in enumerate(STEPS):
    _, light, dark = STAGES[stage]
    cy = ys[i]
    x = CX - BW / 2
    add(f'<rect x="{x}" y="{cy - BH / 2}" width="{BW}" height="{BH}" rx="6" fill="#ffffff" stroke="{dark}" stroke-width="1.3"/>')
    add(f'<circle cx="{x + 22}" cy="{cy}" r="11" fill="{dark}"/>')
    text(x + 22, cy + 4, str(i + 1), 11.5, "#ffffff", "bold")
    text(x + 42, cy - 3, title, 13, INK, "bold", "start")
    text(x + 42, cy + 13, detail, 10.2, MUTED, "normal", "start")
    if i in PARALLEL:
        tag(x + BW - 10, cy - 12, dark)
    if i:
        arrow([(CX, ys[i - 1] + BH / 2), (CX, cy - BH / 2 - 1.5)])

right = CX + BW / 2
left = CX - BW / 2

# design rounds: back from locus selection to primer design
loop_x = PANEL_X - 14              # left of the stage panels, clear of the stage badges
arrow([(left, ys[I_SELECT]), (loop_x, ys[I_SELECT]), (loop_x, ys[I_DESIGN]), (left - 1.5, ys[I_DESIGN])])
mid = (ys[I_SELECT] + ys[I_DESIGN]) / 2
add(f'<rect x="{loop_x - 7}" y="{mid - 64}" width="14" height="128" fill="#ffffff"/>')
text(loop_x, mid, "next round until snp_top_n", 9.5, MUTED, extra=f' dominant-baseline="central" '
     f'transform="rotate(-90 {loop_x} {mid:.1f})"')

# side output of scoring
i_score = 4
document(SIDE_X, ys[i_score], 170, 46, "Locus scores", "SNP summary, plots", DATA_OUT)
arrow([(right, ys[i_score]), (SIDE_X - 86, ys[i_score])])

# optional inputs
i_bg, i_amp = 7, 10
document(SIDE_X, ys[i_bg], 170, 50, "Background genomes", "host, contaminants, community", DATA_IN, dashed=True, stack=True)
arrow([(SIDE_X - 86, ys[i_bg]), (right + 2, ys[i_bg])], dashed=True)
branch_x = right + 30               # background genomes also feed the amplicon screen
arrow([(branch_x, ys[i_bg]), (branch_x, ys[i_amp] - 10), (right + 2, ys[i_amp] - 10)], dashed=True)
database(SIDE_X, ys[i_amp] - 6, 150, 52, "NCBI nt", "remote megablast (default)", DATA_IN)
arrow([(SIDE_X - 76, ys[i_amp] + 10), (right + 2, ys[i_amp] + 10)])

# results (centred on the page, four side by side)
out_y = ys[-1] + 110
RX = W / 2
results = [(RX - 255, "Primer pairs", "sequences, Tm, sites"),
           (RX - 85, "Predicted amplicons", "per genome, FASTA"),
           (RX + 85, "Pooling conflicts", "cross-pair products, CSV"),
           (RX + 255, "Risk taxa", "similar organisms, CSV")]
for x, title, detail in results:
    sx = CX + max(-BW / 2 + 40, min(BW / 2 - 40, (x - CX) * 0.7))
    arrow([(sx, ys[-1] + BH / 2), (x, out_y - 25)])
    document(x, out_y, 162, 46, title, detail, DATA_OUT)

# legend
ly = out_y + 62
lx = PANEL_X
add(f'<rect x="{lx}" y="{ly - 9}" width="16" height="12" rx="1.5" fill="{DATA_IN[0]}" stroke="{DATA_IN[1]}" stroke-width="1.2"/>')
text(lx + 22, ly + 1, "Input data", 10.5, MUTED, anchor="start")
lx += 110
add(f'<rect x="{lx}" y="{ly - 9}" width="16" height="12" rx="1.5" fill="{DATA_OUT[0]}" stroke="{DATA_OUT[1]}" stroke-width="1.2"/>')
text(lx + 22, ly + 1, "Result", 10.5, MUTED, anchor="start")
lx += 90
add(f'<path d="M{lx},{ly - 3} H{lx + 26}" stroke="{ARROW}" stroke-width="1.3" stroke-dasharray="4,3"/>')
text(lx + 32, ly + 1, "Optional input", 10.5, MUTED, anchor="start")
lx += 128
tag(lx + 78, ly - 3, MUTED)
text(lx + 84, ly + 1, "runs in parallel, up to max_cores", 10.5, MUTED, anchor="start")

add("</svg>")
svg[0] = svg[0].replace(f'height="{H}"', f'height="{ly + 20:.0f}"').replace(f'0 0 {W} {H}"', f'0 0 {W} {ly + 20:.0f}"')
svg[3] = svg[3].replace(f'height="{H}"', f'height="{ly + 20:.0f}"')

out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent.parent / "PRISTINE-1.svg"
out.write_text("\n".join(svg) + "\n")
print(f"wrote {out}")
