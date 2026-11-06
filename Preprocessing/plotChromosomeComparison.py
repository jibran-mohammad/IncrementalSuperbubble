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
                       an HPRC chromosome: hprc11, hprcX, ...; HPRC release 2
                       (dataset HPRC2): hprc2-22, hprc2-X, ...              (required)
    -a / --algorithms  algorithm IDs (ALGORITHMS.md); default: every
                       algorithm with data in this mode, in the order of
                       the tables (perEdge: static ... incMinIBSTH1,
                       perBatch: static ... incBSTBatch)
    -o / --output      folder of the PNG        (default: output/plots/comparison)
    -t / --tag         added to the file name, e.g. datastructures for
                       -a incNormalBST incBST incMinIBST
    -s / --steps       only these steps in the bars (default: all), e.g.
                       -s MinI MaxO for the data structures
    -z / --zoom        these algorithms in the zoom panel, on their own axis
                       (default: those below 1/20 of the largest bar);
                       -z none: no zoom panel

More than about six datasets make the bars narrow. Datasets of very
different size (e.g. a small human chromosome next to an HPRC one) are
better plotted separately: the y axis is linear.

=====================================================================
OUTPUT
=====================================================================

    <output>/comparison_<mode>_<datasets>[_<tag>].png
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
from plotcsv import ALL_ALGORITHMS, DISPLAY_NAMES, SHORT_NAMES, STACK_COLORS, ZOOM_RATIO, small_algorithms   # as in the other plots

COMPONENTS = ["UpdateDS", "DFS", "TopologicalOrder", "MinI", "MaxO", "SBCompute"]
MODE_DIR = {"perEdge": "perEdge", "perBatch": "perbatch"}


def dataset_folder(name, mode):
    """short name -> (label, folder with haplotype_build_comparison.csv)"""
    n = name.strip()
    low = n.lower()
    mode_dir = MODE_DIR[mode]
    if low == "ecoli":
        return "Ecoli", REPO / "output" / "Ecoli" / "timeFile" / mode_dir
    m = re.fullmatch(r"hprc2[_-](?:chr)?(\w+)", low)            # HPRC release 2: hprc2-22, hprc2_chrX
    if m:
        c = m.group(1).upper() if m.group(1).isalpha() else m.group(1)
        return f"HPRC v2 chr{c}", REPO / "output" / "HPRC2" / f"chromosome_{c}" / "timeFile" / mode_dir
    m = re.fullmatch(r"hprc[_-]?(?:chr)?(\w+)", low)
    if m:
        c = m.group(1).upper() if m.group(1).isalpha() else m.group(1)
        return f"HPRC chr{c}", REPO / "output" / "HPRC" / f"chromosome_{c}" / "timeFile" / mode_dir
    m = re.fullmatch(r"(?:chr|chromosome_)?(\w+)", low)
    if m:
        c = m.group(1).upper() if m.group(1).isalpha() else m.group(1)
        return f"chr{c}", REPO / "output" / "Human" / f"chromosome_{c}" / "timeFile" / mode_dir
    raise ValueError(f"unknown dataset name: {name}")


def load_steps(folder, algorithms, components, quiet=False):
    """time of every step in components (s), summed over all updates, per algorithm; None if the CSV is missing"""
    csv = folder / "haplotype_build_comparison.csv"
    if not csv.exists():
        print(f"[!] {csv} not found (run parse_to_csv.py), dataset left out")
        return None
    df = pd.read_csv(csv, low_memory=False)
    if "Step" in df.columns:
        df = df[df["Step"].astype(str) != "Base"]
    out = {}
    for a in algorithms:
        # a run that did not finish (stopped at the time limit, or still
        # running: no final PeakMemory line in its time file) has no bar (DNF)
        tf = folder / f"{a}_time.txt"
        if tf.exists():
            with open(tf, "rb") as fh:
                fh.seek(max(0, tf.stat().st_size - 4096))
                if b"PeakMemory" not in fh.read():
                    out[a] = None
                    continue
        cols = [f"{a}{c}Time" for c in components]
        if not any(c in df.columns for c in cols):
            if not quiet:
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
    ap.add_argument("-t", "--tag", help="added to the file name")
    ap.add_argument("-s", "--steps", nargs="+", choices=COMPONENTS, help="only these steps (default: all)")
    ap.add_argument("-z", "--zoom", nargs="+",
                    help="these algorithms in the zoom panel (default: those below 1/ZOOM_RATIO of the largest)")
    args = ap.parse_args()
    algorithms = args.algorithms or list(ALL_ALGORITHMS)
    components = args.steps or COMPONENTS

    labels, steps = [], []
    for name in args.datasets:
        label, folder = dataset_folder(name, args.mode)
        s = load_steps(folder, algorithms, components, quiet=not args.algorithms)
        if s is not None:
            labels.append(label)
            steps.append(s)
    if not args.algorithms:
        # default: only the algorithms that have data in this mode
        algorithms = [a for a in algorithms if any(s[a] is not None for s in steps)]
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
        for ci, comp in enumerate(components):
            for k, a in enumerate(algos):
                vals = np.array([s[a][ci] if s[a] is not None else 0.0 for s in steps])
                ax.bar(groups + offsets[k], vals, width * 0.9, bottom=tops[k], color=STACK_COLORS[comp],
                       edgecolor="black", linewidth=0.5, label=comp if k == 0 else None)
                tops[k] += vals
        pad = 0.015 * tops.max() if tops.max() > 0 else 0
        for k, a in enumerate(algos):
            for i in range(len(labels)):
                name = SHORT_NAMES.get(a, a) if steps[i][a] is not None else f"{SHORT_NAMES.get(a, a)} (DNF)"
                ax.text(groups[i] + offsets[k], tops[k][i] + pad, name, ha="center", va="bottom",
                        fontsize=10, fontweight="bold", rotation=90)
        ax.set_xticks(groups)
        ax.set_xticklabels(labels)
        ax.set_ylim(0, tops.max() * 1.18 if tops.max() > 0 else 1)
        ax.grid(axis="y", linestyle="--", alpha=0.6)

    # algorithms far below the largest one would be invisible: own panel (zoom)
    totals = {a: sum(sum(s[a]) for s in steps if s[a] is not None) for a in algorithms}
    if args.zoom == ["none"]:
        small = []                                      # no zoom panel
    else:
        small = [a for a in args.zoom if a in totals] if args.zoom else small_algorithms(totals)
    width = max(10, 2.2 * len(labels) + 4)
    if small and len(labels) > 6:
        # many datasets: the zoom panel below the main one (side by side the
        # figure gets too wide for the text width)
        fig, (ax, axz) = plt.subplots(2, 1, figsize=(width, 11))
    elif small:
        fig, (ax, axz) = plt.subplots(1, 2, figsize=(width * 1.6, 6), gridspec_kw={"width_ratios": [3, 2]})
    else:
        fig, ax = plt.subplots(figsize=(width, 6))
    draw(ax, algorithms)
    ax.set_ylabel("Time (s), summed over all updates")
    what = "step" if components == COMPONENTS else " and ".join(components)
    ax.set_title(f"Time per {what} and dataset ({args.mode})")
    legend_ax = ax
    if small:
        draw(axz, small)
        axz.set_title(f"Zoom: {', '.join(SHORT_NAMES.get(a, a) for a in small)} (own scale"
                      + (")" if args.zoom else f"; below 1/{ZOOM_RATIO} of the largest)"),
                      fontsize=12)
        legend_ax = axz
    legend_ax.legend(title="Step", loc="upper left", bbox_to_anchor=(1.005, 1.0))
    note = ",  ".join(f"{SHORT_NAMES.get(a, a)} = {DISPLAY_NAMES.get(a, a)}" for a in algorithms)
    fig.text(0.5, 0.005, note, ha="center", va="bottom", fontsize=10)
    plt.tight_layout(rect=(0, 0.04, 1, 1))

    os.makedirs(args.output, exist_ok=True)
    tag = "_".join(l.replace(" ", "") for l in labels)
    if len(tag) > 60:
        # long names break uploads (Overleaf): group and number of datasets instead
        group = "HPRC" if all(l.startswith("HPRC") for l in labels) else \
                ("small" if not any(l.startswith("HPRC") for l in labels) else "mixed")
        tag = f"{group}_{len(labels)}datasets"
    if args.tag:
        tag += f"_{args.tag}"
    path = os.path.join(args.output, f"comparison_{args.mode}_{tag}.png")
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved plot: {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
