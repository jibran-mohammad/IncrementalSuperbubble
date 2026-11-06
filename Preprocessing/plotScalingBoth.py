#!/usr/bin/env python3
"""
plotScalingBoth.py - per edge and per batch in one figure, real (linear) axes:
total time against the number of inserted edges, on all datasets. Two blocks
(per edge, per batch), one column per class (small, medium, large); in a column
the tools in an upper band and our algorithms in a lower band, each with its
own axis starting at 0. Only measured runs: no point = did not finish.

Times as in plotScaling.py: the sum of the timed steps of our algorithms, the
reading of the graph and the search of the tools.

Run from the repository root:
    python3 Preprocessing/plotScalingBoth.py -o output/plots/preview/comparison_both_scaling_linear.png
"""
import argparse
import os
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "Preprocessing"))
from plotScaling import MODE_DIR, edgelist_edges, history_edges, run_time, tool_time   # noqa: E402

CHROMS = ["22", "18", "17", "X", "11"]
CLASSES = [("Small", [f"Human/chromosome_{c}" for c in CHROMS]),
           ("Medium", [f"HPRC/chromosome_{c}" for c in CHROMS]),
           ("Large", ["HPRC2/chromosome_22"])]
ALGOS = {"perEdge": [("Simple", "static"), ("ImpDFS", "incRecompute"), ("ImpDS", "incMinIBST"), ("ImpSB", "incMinIBSTH1")],
         "perBatch": [("Simple", "static"), ("ImpDFS", "incRecomputeBatch"), ("ImpDS", "incBSTBatch"),
                      ("ImpSB", "incBSTBatchH2")]}
STYLE = {"BubbleGun": ("P", "#000000"), "BubbleFinder": ("X", "#7F3C8D"), "Simple": ("s", "#555555"),
         "ImpDFS": ("D", "#E69F00"), "ImpDS": ("^", "#56B4E9"), "ImpSB": ("*", "#D55E00")}
TOOLS = [("BubbleGun", "bubblegun"), ("BubbleFinder", "bubblefinder")]
LIMIT_H = 48.0


def tick(t, _=None):
    """whole numbers with thousands separators, others as they are (2.5, not 2)"""
    return f"{t:,.0f}" if t == int(t) else f"{t:g}"


def name_of(d):
    c = d.split("_")[-1]
    return f"chr{c} v2" if d.startswith("HPRC2") else f"chr{c}"


def load(mode):
    """dataset -> size, and series name -> [(dataset, size, seconds)]"""
    sizes, series = {}, {n: [] for n, _ in ALGOS[mode] + TOOLS}
    for _, ds in CLASSES:
        for d in ds:
            folder = REPO / "output" / d / "timeFile" / MODE_DIR[mode]
            for n, a in ALGOS[mode]:
                r = run_time(folder / f"{a}_time.txt")
                if r:
                    sizes[d] = r[0] if mode == "perEdge" else (history_edges(d) or edgelist_edges(d) or r[0])
                    series[n].append((d, sizes[d], r[1]))
            for n, t in TOOLS:
                r = tool_time(d, t, mode)
                if r and not r[1] and d in sizes:        # estimated: counts as did not finish
                    series[n].append((d, sizes[d], r[0]))
    return sizes, series


def draw(ax, names, series, in_col, unit):
    vals = []
    for n in names:
        pts = sorted([p for p in series[n] if p[0] in in_col], key=lambda p: p[1])
        if not pts:
            continue
        m, c = STYLE[n]
        y = [p[2] / unit for p in pts]
        vals += y
        ax.plot([p[1] for p in pts], y, marker=m, color=c, markersize=15 if m == "*" else 8, markeredgecolor=c,
                linewidth=2 if n == "ImpSB" else 1.3, linestyle="-" if n == "ImpSB" else ":", zorder=2 if n == "ImpSB" else 3)
    ax.grid(True, alpha=0.35, linestyle="--")
    return vals


def main():
    ap = argparse.ArgumentParser(description="Per edge and per batch in one figure, linear axes")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()
    plt.rcParams.update({"font.size": 12})
    fig = plt.figure(figsize=(16, 14))
    outer = fig.add_gridspec(2, 1, hspace=0.3, left=0.06, right=0.985, top=0.93, bottom=0.09)
    for row, (mode, title) in enumerate((("perEdge", "Per edge"), ("perBatch", "Per batch"))):
        sizes, series = load(mode)
        gs = outer[row].subgridspec(2, 3, height_ratios=[1, 2.3], width_ratios=[1, 1.35, 0.4], hspace=0.06, wspace=0.2)
        for col, (cls, ds) in enumerate(CLASSES):
            ds_in = set(ds)
            xs = [sizes[d] for d in ds if d in sizes]
            # perEdge: small in seconds, the others in hours; the tools did not finish there, no upper band
            hours = mode == "perEdge" and col > 0
            unit = 3600.0 if hours else 1.0
            uname = "h" if hours else "s"
            has_tools = any(p[0] in ds_in for n, _ in TOOLS for p in series[n])
            if has_tools:
                top, bottom = fig.add_subplot(gs[0, col]), fig.add_subplot(gs[1, col])
                v = draw(top, [n for n, _ in TOOLS], series, ds_in, unit)
                top.set_ylim(0, max(v) * 1.15)
                top.spines["bottom"].set_visible(False)
                bottom.spines["top"].set_visible(False)
                top.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
                mark = dict(marker=[(-1, -0.6), (1, 0.6)], markersize=9, linestyle="none", color="black",
                            clip_on=False, mew=1.1)
                top.plot([0, 1], [0, 0], transform=top.transAxes, **mark)
                bottom.plot([0, 1], [1, 1], transform=bottom.transAxes, **mark)
                top.yaxis.set_major_formatter(mticker.FuncFormatter(tick))
                axes = [top, bottom]
            else:
                bottom = fig.add_subplot(gs[:, col])
                axes = [bottom]
            v = draw(bottom, [n for n, _ in ALGOS[mode]], series, ds_in, unit)
            hi = max(v) * 1.15
            if hours:
                hi = LIMIT_H * 1.15                       # the same axis for medium and large, with the limit
                bottom.axhline(LIMIT_H, color="black", linestyle="--", linewidth=2.2)
                bottom.text(0.98, LIMIT_H + 0.8, "time limit (2 days)", transform=bottom.get_yaxis_transform(),
                            ha="right", fontsize=11, fontweight="bold")
            bottom.set_ylim(0, hi)
            bottom.yaxis.set_major_formatter(mticker.FuncFormatter(tick))
            for ax in axes:
                if len(xs) == 1:
                    ax.set_xlim(xs[0] * 0.8, xs[0] * 1.2)
                else:
                    span = max(xs) - min(xs)
                    ax.set_xlim(min(xs) - 0.07 * span, max(xs) + 0.07 * span)
            # one tick per dataset: its name and size, so no name covers a point
            on_axis = sorted((sizes[d], d) for d in ds if d in sizes)
            bottom.set_xticks([x for x, _ in on_axis])
            bottom.set_xticklabels([f"{name_of(d)}\n{x:,}" for x, d in on_axis], rotation=90, fontsize=9, linespacing=1.0)
            for ax in axes[:-1]:
                ax.set_xticks([x for x, _ in on_axis])
            bottom.set_xlabel("inserted edges")
            bottom.set_ylabel(f"total time ({uname})")
            axes[0].set_title(f"{title}: {cls}", fontsize=13)
    handles = [plt.Line2D([], [], marker=STYLE[n][0], color=STYLE[n][1], markersize=15 if STYLE[n][0] == "*" else 8,
                          linestyle="-" if n == "ImpSB" else ":", label=n)
               for n in ["BubbleGun", "BubbleFinder", "Simple", "ImpDFS", "ImpDS", "ImpSB"]]
    fig.legend(handles=handles, loc="upper center", ncol=len(handles), frameon=False, fontsize=12)
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    fig.savefig(args.output, dpi=180)
    print(f"Saved plot: {args.output}")


if __name__ == "__main__":
    sys.exit(main())
