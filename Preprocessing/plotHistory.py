#!/usr/bin/env python3
"""
plotHistory.py - figures from the superbubble history (the summary written
by Preprocessing/superbubbleHistory.py, one row per genome).

=====================================================================
HOW TO RUN
=====================================================================

    # figure 1 (growth curve) and figure 2 (new sites per genome by
    # population) for one or more histories
    python3 Preprocessing/plotHistory.py \\
        --summary output/HPRC/chromosome_22/history/chr22_summary.tsv \\
                  output/HPRC2/chromosome_22/history/chr22_summary.tsv \\
        --label "chr22 release 1" "chr22 release 2" \\
        --out output/plots/history/hprc_chr22

Options:
    --summary   <out>_summary.tsv of superbubbleHistory.py (with
                --metadata for figure 2), one or more
    --label     a label per summary (default: the file name)
    --group     metadata column for the colours of figure 2
                (default: Superpopulation)
    --lifespans, --alleles
                for figure 3: <out>_lifespans.tsv (superbubbleHistory.py
                with --removed-links) and <out>_alleles.tsv
                (superbubbleAlleles.py), one each per summary, in the same
                order
    --out       output prefix

=====================================================================
OUTPUT
=====================================================================

    <out>_growth.png      figure 1: superbubbles alive, cumulative born
                          (new sites: containing a vertex the genome added),
                          cumulative promoted (an existing region became a
                          superbubble) and cumulative died, against the
                          number of genomes added (one panel per summary)
    <out>_new_sites.png   figure 2: new sites (born) per genome in the order
                          the genomes were added, coloured by --group, and
                          their distribution per group (one row per summary
                          with metadata)
    <out>_var_types.png   figure 3 (with --lifespans and --alleles): variant
                          type of the final superbubbles, untouched vs
                          touched by the DAG normalization (dag_affected),
                          as shares of each group, with the counts
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.patches import Patch
import pandas as pd

GROUP_COLORS = {"AFR": "#D55E00", "AMR": "#009E73", "EAS": "#0072B2", "EUR": "#CC79A7",
                "SAS": "#E69F00", "-": "#999999"}


def main():
    ap = argparse.ArgumentParser(description="Figures from the superbubble history")
    ap.add_argument("--summary", nargs="+", required=True)
    ap.add_argument("--label", nargs="+")
    ap.add_argument("--group", default="Superpopulation")
    ap.add_argument("--lifespans", nargs="+")
    ap.add_argument("--alleles", nargs="+")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    labels = args.label or [os.path.basename(p).replace("_summary.tsv", "") for p in args.summary]
    frames = [pd.read_csv(p, sep="\t") for p in args.summary]
    if os.path.dirname(args.out):
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
    plt.rcParams.update({"font.size": 13, "axes.labelsize": 13, "legend.fontsize": 11})

    # ---- figure 1: growth curve
    fig, axes = plt.subplots(1, len(frames), figsize=(9.5 * len(frames), 5), squeeze=False)
    for ax, df, lab in zip(axes[0], frames, labels):
        g = df[df["step"] > 0]
        base = int(df.loc[df["step"] == 0, "alive"].iloc[0]) if (df["step"] == 0).any() else 0
        ax.plot(g["step"], g["alive"], color="black", linewidth=2, label="alive")
        ax.plot(g["step"], base + g["cum_born"], color="#0072B2", label="base + cumulative born")
        ax.plot(g["step"], base + g["cum_born"] + g["cum_promoted"], color="#009E73", linestyle="--",
                label="... + cumulative promoted")
        ax.plot(g["step"], g["cum_died"], color="#D55E00", label="cumulative died")
        ax.axhline(base, color="grey", linestyle=":", linewidth=1, label=f"base graph ({base:,})")
        ax.set_xlabel("genomes added")
        ax.set_ylabel("superbubbles")
        ax.set_title(lab)
        ax.grid(alpha=0.4, linestyle="--")
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
        ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))
    plt.tight_layout()
    plt.savefig(args.out + "_growth.png", dpi=200)
    plt.close(fig)
    print(f"Saved plot: {args.out}_growth.png")

    # ---- figure 3: variant types, untouched vs touched by the DAG normalization
    if args.lifespans and args.alleles:
        plot_var_types(args, labels)

    # ---- figure 2: new sites per genome by population
    with_meta = [(df, lab) for df, lab in zip(frames, labels) if args.group in df.columns]
    if not with_meta:
        print(f"[!] no column {args.group} (run superbubbleHistory.py with --metadata): figure 2 skipped")
        return
    fig, axes = plt.subplots(len(with_meta), 2, figsize=(16, 4.5 * len(with_meta)), squeeze=False,
                             gridspec_kw={"width_ratios": [3, 1]})
    for (ax, axb), (df, lab) in zip(axes, with_meta):
        g = df[df["step"] > 0].copy()
        g[args.group] = g[args.group].fillna("-").astype(str)
        colors = [GROUP_COLORS.get(v, "#555555") for v in g[args.group]]
        ax.bar(g["step"], g["born"], color=colors, width=0.9)
        ax.set_xlabel("genome (in the order added)")
        ax.set_ylabel("new sites (born)")
        ax.set_title(f"{lab}: new superbubble sites per genome")
        ax.grid(axis="y", alpha=0.4, linestyle="--")
        groups = [k for k in GROUP_COLORS if k in set(g[args.group])] + \
                 sorted(set(g[args.group]) - set(GROUP_COLORS))
        handles = [Patch(color=GROUP_COLORS.get(k, "#555555"), label=f"{k} ({(g[args.group] == k).sum()})")
                   for k in groups]
        ax.legend(handles=handles, title=f"{args.group} (genomes)", ncol=len(groups), loc="upper right")
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
        # distribution per group, without the first genomes (they add the most sites to any graph)
        late = g[g["step"] > max(3, len(g) // 10)]
        data = [late.loc[late[args.group] == k, "born"] for k in groups]
        bp = axb.boxplot(data, patch_artist=True, showfliers=False)
        for patch, k in zip(bp["boxes"], groups):
            patch.set_facecolor(GROUP_COLORS.get(k, "#555555"))
        for i, d in enumerate(data, 1):
            axb.scatter([i] * len(d), d, color="black", s=8, zorder=3)
        axb.set_xticks(range(1, len(groups) + 1))
        axb.set_xticklabels(groups)
        axb.set_title(f"per {args.group} (genomes {int(late['step'].min())}+)", fontsize=11)
        axb.grid(axis="y", alpha=0.4, linestyle="--")
    plt.tight_layout()
    plt.savefig(args.out + "_new_sites.png", dpi=200)
    plt.close(fig)
    print(f"Saved plot: {args.out}_new_sites.png")


VAR_ORDER = ["SNP", "MNP", "insertion", "deletion", "complex", "mixed", "nested", "nonref", "none"]
VAR_COLORS = {"SNP": "#0072B2", "MNP": "#56B4E9", "insertion": "#009E73", "deletion": "#D55E00",
              "complex": "#CC79A7", "mixed": "#E69F00", "nested": "#999999", "nonref": "#000000",
              "none": "#DDDDDD"}


def plot_var_types(args, labels):
    """figure 3: one bar per dataset and dag_affected group, stacked by variant type (shares)"""
    bars = []                                      # (label, counts per type, total)
    for lab, lf, af in zip(labels, args.lifespans, args.alleles):
        life = pd.read_csv(lf, sep="\t", usecols=["entrance", "exit", "died_step", "dag_affected"],
                           dtype=str, keep_default_na=False)          # "NA" = alive at the end, a value
        life = life[life["died_step"] == "NA"]
        al = pd.read_csv(af, sep="\t", dtype=str, keep_default_na=False)
        df = al.merge(life[["entrance", "exit", "dag_affected"]], on=["entrance", "exit"], how="left")
        for flag, name in (("0", "untouched"), ("1", "touched by DAG")):
            sub = df[df["dag_affected"] == flag]
            if len(sub):
                bars.append((f"{lab}\n{name}", sub["var_type"].value_counts(), len(sub)))
    if not bars:
        print("[!] no dag_affected values (superbubbleHistory.py --removed-links): figure 3 skipped")
        return
    fig, ax = plt.subplots(figsize=(max(8, 1.6 * len(bars) + 4), 6))
    for i, (lab, counts, total) in enumerate(bars):
        bottom = 0.0
        for t in VAR_ORDER:
            share = 100.0 * counts.get(t, 0) / total
            ax.bar(i, share, bottom=bottom, color=VAR_COLORS[t], edgecolor="black", linewidth=0.4,
                   label=t if i == 0 else None)
            bottom += share
        ax.text(i, 101, f"{total:,}", ha="center", va="bottom", fontsize=10)
    ax.set_xticks(range(len(bars)))
    ax.set_xticklabels([b[0] for b in bars], fontsize=10)
    ax.set_ylabel("share of the final superbubbles (%)")
    ax.set_ylim(0, 108)
    ax.set_title("Variant type of the final superbubbles (number above each bar)")
    ax.legend(title="var_type", loc="upper left", bbox_to_anchor=(1.01, 1.0))
    plt.tight_layout()
    plt.savefig(args.out + "_var_types.png", dpi=200)
    plt.close(fig)
    print(f"Saved plot: {args.out}_var_types.png")


if __name__ == "__main__":
    main()
