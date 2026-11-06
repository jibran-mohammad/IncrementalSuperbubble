#!/usr/bin/env python3
"""
plotGroups.py - one row of panels, one panel per group of datasets (e.g. the
smaller human chromosomes, the other human chromosomes, HPRC chr22), in each
panel one group of bars per dataset and one bar per algorithm, each bar
stacked from the time of each step (UpdateDS ... SBCompute) summed over all
updates. Each panel has its own linear axis; with --log all panels share one
logarithmic axis (comparable across the groups).

=====================================================================
HOW TO RUN
=====================================================================

    python3 Preprocessing/plotGroups.py -m perEdge \\
        --group "chr13 / 14 / 15 / 21" 13 14 15 21 \\
        --group "chr11 / 17 / 18 / 22 / X" 11 17 18 22 X \\
        --group "HPRC chr22" hprc22 \\
        -a static incRecompute incMinIBST incMinIBSTH1 \\
        -o output/plots/comparison/groups/comparison_perEdge_groups.png [--log]

Options:
    -m / --mode        perEdge | perBatch
    --group TITLE D..  a panel: its title, then its datasets (short names as
                       in plotChromosomeComparison.py: 13, chrX, hprc22, ...)
    -a / --algorithms  algorithm IDs (one bar each, in this order)
    --log              one shared logarithmic axis instead of one linear axis
                       per panel
    -o / --output      the PNG
"""
import argparse
import os
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "Preprocessing"))
from plotcsv import SHORT_NAMES, STACK_COLORS                                       # noqa: E402
from plotChromosomeComparison import COMPONENTS, dataset_folder, load_steps          # noqa: E402


def main():
    ap = argparse.ArgumentParser(description="Groups of datasets in one row of panels")
    ap.add_argument("-m", "--mode", required=True, choices=["perEdge", "perBatch"])
    ap.add_argument("--group", nargs="+", action="append", required=True, metavar="TITLE DATASET")
    ap.add_argument("-a", "--algorithms", nargs="+", required=True)
    ap.add_argument("--log", action="store_true")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    groups = []                                     # (title, [(label, steps)])
    for g in args.group:
        title, names = g[0], g[1:]
        data = []
        for n in names:
            label, folder = dataset_folder(n, args.mode)
            s = load_steps(folder, args.algorithms, COMPONENTS, quiet=True)
            if s is not None:
                data.append((label, s))
        groups.append((title, data))

    plt.rcParams.update({"font.size": 13, "axes.labelsize": 13, "xtick.labelsize": 12,
                         "ytick.labelsize": 11, "legend.fontsize": 11})
    widths = [max(1.2, len(d)) for _, d in groups]
    fig, axes = plt.subplots(1, len(groups), figsize=(3.2 * sum(widths) / 2 + 3, 5.2), sharey=args.log,
                             gridspec_kw={"width_ratios": widths})
    axes = np.atleast_1d(axes)
    algos = args.algorithms
    n = len(algos)
    width = 0.8 / n
    offsets = (np.arange(n) - (n - 1) / 2) * width
    for ax, (title, data) in zip(axes, groups):
        pos = np.arange(len(data))
        tops = np.zeros((n, len(data)))
        for ci, comp in enumerate(COMPONENTS):
            for k, a in enumerate(algos):
                vals = np.array([s[a][ci] if s[a] is not None else 0.0 for _, s in data])
                ax.bar(pos + offsets[k], vals, width * 0.9, bottom=tops[k], color=STACK_COLORS[comp],
                       edgecolor="black", linewidth=0.4, label=comp)
                tops[k] += vals
        for k, a in enumerate(algos):
            for i, (_, s) in enumerate(data):
                name = SHORT_NAMES.get(a, a) if s[a] is not None else f"{SHORT_NAMES.get(a, a)} (DNF)"
                y = tops[k][i] if tops[k][i] > 0 else (1 if args.log else 0)
                ax.text(pos[i] + offsets[k], y * (1.08 if args.log else 1) + (0 if args.log else 0.015 * tops.max()),
                        name, ha="center", va="bottom", fontsize=8.5, fontweight="bold", rotation=90)
        ax.set_xticks(pos)
        ax.set_xticklabels([l.replace("HPRC ", "").replace("v2 chr22", "chr22 v2") for l, _ in data])
        ax.set_title(title, fontsize=13)
        ax.grid(axis="y", linestyle="--", alpha=0.5)
        if not args.log:
            ax.set_ylim(0, tops.max() * 1.25 if tops.max() > 0 else 1)
    if args.log:
        axes[0].set_yscale("log")
        allv = [t for ax in axes for b in ax.patches for t in [b.get_y() + b.get_height()] if t > 0]
        axes[0].set_ylim(min(allv) / 3, max(allv) * 6)
    axes[0].set_ylabel("Time (s), summed over all updates" + (" (log scale)" if args.log else ""))
    handles, labels = axes[0].get_legend_handles_labels()
    seen = {}
    for h, l in zip(handles, labels):
        seen.setdefault(l, h)
    fig.legend(list(seen.values()), list(seen.keys()), loc="lower center", ncol=len(seen), frameon=False)
    names = ",  ".join(f"{SHORT_NAMES.get(a, a)} = {a}" for a in algos)
    plt.tight_layout(rect=(0, 0.08, 1, 1))
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    fig.savefig(args.output, dpi=250)
    print(f"Saved plot: {args.output}")


if __name__ == "__main__":
    sys.exit(main())
