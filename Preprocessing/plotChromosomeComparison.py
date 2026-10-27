#!/usr/bin/env python3
r"""
plotChromosomeComparison.py

One bar chart comparing algorithms across datasets: for every dataset one
group of bars, one bar per algorithm, each bar stacked from the time of
each step (UpdateDS, DFS, TopologicalOrder, MinI, MaxO, SBCompute) summed
over all updates (the base graph is left out). Works for perEdge and
perBatch; the datasets are chosen by short names.

=====================================================================
HOW TO RUN
=====================================================================

Run from the repository root, with the virtual environment active, after
Preprocessing/parse_to_csv.py for every dataset in the plot.

    # perEdge, E. coli and four human chromosomes
    python3 Preprocessing/plotChromosomeComparison.py -m perEdge -d ecoli 11 13 22 X

    # perBatch, HPRC chromosomes
    python3 Preprocessing/plotChromosomeComparison.py -m perBatch -d hprc11 hprc18 hprcX

    # other algorithms
    python3 Preprocessing/plotChromosomeComparison.py -m perEdge -d 11 22 X -a static incBST incMinIBSTH1

Options:

    -m / --mode        perEdge | perBatch                                   (required)
    -d / --datasets    ecoli; a human chromosome: 11, chr11, X, chrX, ...;
                       an HPRC chromosome: hprc11, hprcX, ...                (required)
    -a / --algorithms  algorithm IDs (ALGORITHMS.md); default perEdge:
                       static incRecompute incMinIBSTH1, perBatch: static
                       incMinIBSTH1 incBSTBatch
    -o / --output      folder of the PNG        (default: output/plots/comparison)

More than about six datasets make the bars narrow. Datasets of very
different size (e.g. a small human chromosome next to an HPRC one) are
better plotted separately: the y axis is linear.

=====================================================================
OUTPUT
=====================================================================

    <output>/comparison_<mode>_<datasets>.png
"""
import argparse
import os
import pathlib
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "Preprocessing"))
from plotcsv import DISPLAY_NAMES, SHORT_NAMES, STACK_COLORS, ZOOM_RATIO, small_algorithms   # as in the other plots

COMPONENTS = ["UpdateDS", "DFS", "TopologicalOrder", "MinI", "MaxO", "SBCompute"]
DEFAULT_ALGORITHMS = {"perEdge": ["static", "incRecompute", "incMinIBSTH1"],
                      "perBatch": ["static", "incMinIBSTH1", "incBSTBatch"]}
MODE_DIR = {"perEdge": "perEdge", "perBatch": "perbatch"}


def dataset_folder(name, mode):
    """short name -> (label, folder with haplotype_build_comparison.csv)"""
    n = name.strip()
    low = n.lower()
    mode_dir = MODE_DIR[mode]
    if low == "ecoli":
        return "Ecoli", REPO / "output" / "Ecoli" / "timeFile" / mode_dir
    m = re.fullmatch(r"hprc[_-]?(?:chr)?(\w+)", low)
    if m:
        c = m.group(1).upper() if m.group(1).isalpha() else m.group(1)
        return f"HPRC chr{c}", REPO / "output" / "HPRC" / f"chromosome_{c}" / "timeFile" / mode_dir
    m = re.fullmatch(r"(?:chr|chromosome_)?(\w+)", low)
    if m:
        c = m.group(1).upper() if m.group(1).isalpha() else m.group(1)
        return f"chr{c}", REPO / "output" / "Human" / f"chromosome_{c}" / "timeFile" / mode_dir
    raise ValueError(f"unknown dataset name: {name}")


def load_steps(folder, algorithms):
    """time of every step (s), summed over all updates, per algorithm; None if the CSV is missing"""
    csv = folder / "haplotype_build_comparison.csv"
    if not csv.exists():
        print(f"[!] {csv} not found (run parse_to_csv.py), dataset left out")
        return None
    df = pd.read_csv(csv, low_memory=False)
    if "Step" in df.columns:
        df = df[df["Step"].astype(str) != "Base"]
    out = {}
    for a in algorithms:
        cols = [f"{a}{c}Time" for c in COMPONENTS]
        if not any(c in df.columns for c in cols):
            print(f"[!] {csv}: no data of {a}")
            out[a] = None
            continue
        out[a] = [float(pd.to_numeric(df[c], errors="coerce").fillna(0).sum()) if c in df.columns else 0.0
                  for c in cols]
    return out


def main():
    ap = argparse.ArgumentParser(description="Stacked step times of several algorithms on several datasets")
    ap.add_argument("-m", "--mode", required=True, choices=["perEdge", "perBatch"])
    ap.add_argument("-d", "--datasets", nargs="+", required=True)
    ap.add_argument("-a", "--algorithms", nargs="+")
    ap.add_argument("-o", "--output", default=str(REPO / "output" / "plots" / "comparison"))
    args = ap.parse_args()
    algorithms = args.algorithms or DEFAULT_ALGORITHMS[args.mode]

    labels, steps = [], []
    for name in args.datasets:
        label, folder = dataset_folder(name, args.mode)
        s = load_steps(folder, algorithms)
        if s is not None:
            labels.append(label)
            steps.append(s)
    if not labels:
        print("No data found.")
        return 1

    plt.rcParams.update({"font.size": 14, "axes.labelsize": 14, "xtick.labelsize": 13,
                         "ytick.labelsize": 13, "legend.fontsize": 12, "legend.title_fontsize": 13})

    def draw(ax, algos):
        """one group per dataset, one stacked bar per algorithm, its name above it"""
        n = len(algos)
        width = 0.8 / n
        groups = np.arange(len(labels))
        offsets = (np.arange(n) - (n - 1) / 2) * width
        tops = np.zeros((n, len(labels)))
        for ci, comp in enumerate(COMPONENTS):
            for k, a in enumerate(algos):
                vals = np.array([s[a][ci] if s[a] is not None else 0.0 for s in steps])
                ax.bar(groups + offsets[k], vals, width * 0.9, bottom=tops[k], color=STACK_COLORS[comp],
                       edgecolor="black", linewidth=0.5, label=comp if k == 0 else None)
                tops[k] += vals
        pad = 0.015 * tops.max() if tops.max() > 0 else 0
        for k, a in enumerate(algos):
            for i in range(len(labels)):
                name = SHORT_NAMES.get(a, a) if steps[i][a] is not None else f"{SHORT_NAMES.get(a, a)} (-)"
                ax.text(groups[i] + offsets[k], tops[k][i] + pad, name, ha="center", va="bottom",
                        fontsize=10, fontweight="bold", rotation=90)
        ax.set_xticks(groups)
        ax.set_xticklabels(labels)
        ax.set_ylim(0, tops.max() * 1.18 if tops.max() > 0 else 1)
        ax.grid(axis="y", linestyle="--", alpha=0.6)

    # algorithms far below the largest one would be invisible: own panel (zoom)
    totals = {a: sum(sum(s[a]) for s in steps if s[a] is not None) for a in algorithms}
    small = small_algorithms(totals)
    width = max(10, 2.2 * len(labels) + 4)
    if small:
        fig, (ax, axz) = plt.subplots(1, 2, figsize=(width * 1.6, 6), gridspec_kw={"width_ratios": [3, 2]})
    else:
        fig, ax = plt.subplots(figsize=(width, 6))
    draw(ax, algorithms)
    ax.set_ylabel("Time (s), summed over all updates")
    ax.set_title(f"Time per step and dataset ({args.mode})")
    legend_ax = ax
    if small:
        draw(axz, small)
        axz.set_title(f"Zoom: {', '.join(SHORT_NAMES.get(a, a) for a in small)} (own scale; below 1/{ZOOM_RATIO} of the largest)",
                      fontsize=12)
        legend_ax = axz
    legend_ax.legend(title="Step", loc="upper left", bbox_to_anchor=(1.005, 1.0))
    note = ",  ".join(f"{SHORT_NAMES.get(a, a)} = {DISPLAY_NAMES.get(a, a)}" for a in algorithms)
    fig.text(0.5, 0.005, note, ha="center", va="bottom", fontsize=10)
    plt.tight_layout(rect=(0, 0.04, 1, 1))

    os.makedirs(args.output, exist_ok=True)
    tag = "_".join(l.replace(" ", "") for l in labels)
    path = os.path.join(args.output, f"comparison_{args.mode}_{tag}.png")
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved plot: {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
