#!/usr/bin/env python3

"""
plotChromosomeComparison.py

Generate a bar plot comparing algorithm runtime across chromosomes.

The script expects directories that contain the file:

    haplotype_build_comparison.csv

Example directory:

/home/.../chromosome_11/timeFile/perEdge/
    haplotype_build_comparison.csv


Example usage:

python plotChromosomeComparison.py -p path/chromosome_1/timeFile/perEdge path/chromosome_2/timeFile/perEdge \
   path/chromosome_11/timeFile/perEdge \
-m perEdge \
-o plots


Maximum number of chromosomes allowed: 5
"""


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os
import re


ALGORITHMS = ["static", "incBSTBatch", "incMinIBST"]

DISPLAY_NAMES = {
    "static": "Static",
    "incBSTBatch": "Inc-BST-Batch",
    "incMinIBST": "FINAL"
}

def plot_component_comparison(paths, output_dir, mode):

    plt.rcParams.update({
        "font.size": 14,
        "axes.titlesize": 16,
        "axes.labelsize": 14,
        "xtick.labelsize": 13,
        "ytick.labelsize": 13,
        "legend.fontsize": 13,
        "legend.title_fontsize": 14,
    })

    COMPONENTS = [
        "Insert",
        "DFS",
        "TopologicalOrder",
        "MinI",
        "MaxO",
        "SBCompute"
    ]

    labels = []
    totals = {}

    # structure: totals[algorithm][component]
    for alg in ALGORITHMS:
        totals[alg] = {comp: [] for comp in COMPONENTS}

    for path in paths:

        label = extract_chromosome_label(path)
        labels.append(label)

        df = load_csv(path)

        for alg in ALGORITHMS:

            for comp in COMPONENTS:

                col = f"{alg}{comp}Time"

                value = df[col].sum()

                totals[alg][comp].append(value)

    # plotting
    bar_width = 0.17
    step_indices = np.arange(len(labels))

    x_positions = {
        "static":   step_indices - bar_width,
        "incBSTBatch": step_indices ,
        "incMinIBST": step_indices +  bar_width,
    }

    fig, ax = plt.subplots(figsize=(20,6))

    version_bottoms = {
        alg: np.zeros(len(labels))
        for alg in ALGORITHMS
    }

    STACK_COLORS = {
        "Insert": "#4C72B0",
        "DFS": "#55A868",
        "TopologicalOrder": "#DD8452",
        "MinI": "#C44E52",
        "MaxO": "#8172B2",
        "SBCompute": "#CCB974"
    }

    for comp in COMPONENTS:

        color = STACK_COLORS[comp]

        for alg in ALGORITHMS:

            data = totals[alg][comp]

            label = comp if alg == "static" else None

            ax.bar(
                x_positions[alg],
                data,
                bar_width,
                bottom=version_bottoms[alg],
                color=color,
                edgecolor="black",
                linewidth=0.5,
                label=label
            )

            version_bottoms[alg] += np.array(data)

    # add static / incRecompute / incBST labels
    for i in range(len(labels)):
        for alg in ALGORITHMS:

            x = x_positions[alg][i]
            y = version_bottoms[alg][i]

            ax.text(
                x,   
                y + 0.02 * max(max(v) for v in version_bottoms.values()),
                DISPLAY_NAMES.get(alg, alg),
                ha="center",
                fontsize=11,
                fontweight="bold",
                rotation=90,
            )

    ax.set_xticks(step_indices)
    ax.set_xticklabels(labels)

    ax.set_ylabel("Time (seconds)")
    #ax.set_title(f"Chromosome Runtime Breakdown ({mode})")

    ax.legend(title="Component")

    ax.grid(axis="y", linestyle="--", alpha=0.6)

    plt.tight_layout()

    filename = build_filename(mode, labels)

    full_path = os.path.join(output_dir, filename)

    plt.savefig(full_path, dpi=300)

    plt.close()

    print(f"Saved plot: {full_path}")


def extract_chromosome_label(path):

    # Ecoli detection
    if "Ecoli" in path or "ecoli" in path:
        return "Ecoli"

    # chromosome numbers
    match = re.search(r"chromosome_(\d+)", path)
    if match:
        return f"chr{match.group(1)}"

    # chromosome X
    match = re.search(r"chromosome_X", path)
    if match:
        return "chrX"

    # fallback
    return os.path.basename(os.path.normpath(path))


def load_csv(directory):

    csv_path = os.path.join(directory, "haplotype_build_comparison.csv")

    if not os.path.exists(csv_path):
        print(f"Warning: {csv_path} not found")
        return None

    print(f"Reading {csv_path}")

    df = pd.read_csv(csv_path)

    if "Step" in df.columns:
        df["Step"] = df["Step"].replace({"Base": 0})
        df["Step"] = pd.to_numeric(df["Step"], errors="coerce").fillna(0)

    df = df[df["Step"] != 0]

    for alg in ALGORITHMS:
        col = f"{alg}TotalTime"
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)

    return df


def compute_totals(paths):

    totals = {alg: [] for alg in ALGORITHMS}
    labels = []

    for path in paths:

        df = load_csv(path)

        if df is None:
            continue

        label = extract_chromosome_label(path)

        labels.append(label)

        for alg in ALGORITHMS:
            totals[alg].append(df[f"{alg}TotalTime"].sum())

    return totals, labels


def build_filename(mode, labels):

    chrom_string = "_".join(labels)

    alg_string = "_".join([a.lower() for a in ALGORITHMS])

    filename = f"comparison_{mode}_{chrom_string}_{alg_string}.png"

    return filename




def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-p",
        "--paths",
        nargs="+",
        required=True,
        help="Directories containing haplotype_build_comparison.csv"
    )

    parser.add_argument(
        "-m",
        "--mode",
        choices=["perEdge", "perbatch"],
        required=True,
        help="Experiment mode"
    )

    parser.add_argument(
        "-o",
        "--output",
        default="plots",
        help="Output directory"
    )

    args = parser.parse_args()

    totals, labels = compute_totals(args.paths)

    if len(labels) == 0:
        print("No valid data found")
        return

    filename = build_filename(args.mode, labels)

    plot_component_comparison(
        args.paths,
        args.output,
        args.mode
    )


if __name__ == "__main__":
    main()