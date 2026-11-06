#!/usr/bin/env python3
"""
plotScaling.py - total time of the algorithms against the size of the
dataset (number of inserted edges), one point per finished run, log scales,
with the time limit as a line. perEdge: shows which algorithms scale to the
large graphs (H1 on all of them) and where the others stop (DNF: no point);
perBatch: how the batch algorithms and the tools grow with the dataset.

=====================================================================
HOW TO RUN
=====================================================================

    python3 Preprocessing/plotScaling.py \\
        --datasets HPRC/chromosome_11 HPRC/chromosome_17 HPRC/chromosome_18 \\
                   HPRC/chromosome_22 HPRC/chromosome_X HPRC2/chromosome_22 \\
        --out output/plots/comparison/medium/comparison_perEdge_HPRC_scaling.png

    # all datasets in one plot, with broken axes (the two plots of the paper)
    python3 Preprocessing/plotScaling.py --broken --tools bubblegun bubblefinder \\
        --datasets Human/chromosome_22 Human/chromosome_18 Human/chromosome_17 Human/chromosome_11 Human/chromosome_X \\
                   HPRC/chromosome_22 HPRC/chromosome_18 HPRC/chromosome_17 HPRC/chromosome_X HPRC/chromosome_11 \\
                   HPRC2/chromosome_22 \\
        --out output/plots/comparison/groups/comparison_perEdge_scaling.png

    # perBatch: the batch algorithms and the tools (measured, one update per haplotype)
    python3 Preprocessing/plotScaling.py --mode perBatch \\
        --datasets HPRC/chromosome_22 HPRC/chromosome_18 HPRC/chromosome_17 \\
                   HPRC/chromosome_X HPRC/chromosome_11 HPRC2/chromosome_22 \\
        --algorithms static incRecomputeBatch incBSTBatch incBSTBatchH1 \\
        --out output/plots/comparison/medium/comparison_perBatch_HPRC_scaling.png

Options:
    --mode       perEdge (default; time in hours) | perBatch (seconds; the
                 size is the number of inserted edges over all haplotypes,
                 from the history summary of the dataset)
    --split      two panels: our algorithms on a linear axis, the tools (with the best algorithm as
                 the reference) on a logarithmic axis, for when the tools are far above the algorithms
    --broken     one plot with one column per class of datasets (perEdge: small, medium + large;
                 perBatch: small, medium, large), the tools in an upper band and our algorithms in a
                 lower band, each part with its own log range, so the differences stay visible. Only
                 measured runs: estimates count as did not finish
    --datasets   folders under output/ with timeFile/<mode>/<ID>_time.txt
    --algorithms algorithm IDs (default: static incRecompute incNormalBST
                 incBST incMinIBST incMinIBSTH1)
    --limit-days the time limit drawn as a line (default 2)
    --out        the PNG

A run counts only if it finished (its time file ends with the PeakMemory
line); its time is the sum of the timed steps of every update (Total_Time),
as in the time tables.
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
from plotcsv import SHORT_NAMES                                      # noqa: E402

TOOL_NAMES = {"bubblegun": "BubbleGun", "bubblefinder": "BubbleFinder"}
MARKERS = {"bubblegun": ("P", "#000000"), "bubblefinder": ("X", "#7F3C8D"),
           "static": ("s", "#555555"), "incRecompute": ("D", "#E69F00"), "incNormalBST": ("v", "#CC79A7"),
           "incBST": ("^", "#56B4E9"), "incMinIBST": ("o", "#009E73"), "incMinIBSTH1": ("*", "#D55E00"),
           "incRecomputeBatch": ("D", "#E69F00"), "incBSTBatch": ("^", "#56B4E9"), "incBSTBatchH1": ("o", "#009E73"), "incBSTBatchH2": ("*", "#D55E00")}
MODE_DIR = {"perEdge": "perEdge", "perBatch": "perbatch"}


def run_time(path):
    """(inserted edges, total time in s) of a finished run, else None"""
    if not path.exists():
        return None
    with open(path, "rb") as fh:
        fh.seek(max(0, path.stat().st_size - 4096))
        if b"PeakMemory" not in fh.read():
            return None
    edges, total = 0, 0.0
    with open(path) as f:
        for line in f:
            p = line.split()
            # update lines: "<number> <genome> <UpdateDS> ... <Total_Time> ..."
            if len(p) >= 10 and p[0].isdigit() and not p[1].replace(".", "").isdigit():
                edges += 1
                total += float(p[8])
    return edges, total


def label(dataset):
    rel, chrom = dataset.split("/")
    c = chrom.replace("chromosome_", "chr")
    return f"{c} v2" if rel == "HPRC2" else c


def history_edges(dataset):
    """inserted edges over all haplotypes (the base graph excluded) from the
    history summary of the dataset (new_edges of the steps after step 0), else None"""
    for h in sorted((REPO / "output" / dataset / "history").glob("*_history_summary.tsv")):
        with open(h) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            if "new_edges" not in head:
                continue
            ce, cs = head.index("new_edges"), head.index("step")
            return sum(int(l.split("\t")[ce]) for l in fh if l.split("\t")[cs] != "0")
    return None


def edgelist_edges(dataset):
    """inserted edges over all haplotypes (the base graph excluded) from the edge list, else None"""
    f = REPO / "input" / dataset / "Graphs" / "numeric_edgelist.txt"
    if not f.exists():
        return None
    n, started = 0, False
    with open(f) as fh:
        for line in fh:
            if line.startswith("# HAPLOTYPE"):
                started = True
            elif started and not line.startswith("#") and len(line.split()) == 2:
                n += 1
    return n


def tool_time(dataset, tool, mode="perEdge"):
    """(total time in s, estimated) of an external tool, from
    output/<dataset>/<tool>/<mode>/<tool>_time.txt (TotalTime; estimated if
    it was timed on a sample of the updates), else None"""
    f = REPO / "output" / dataset / tool / MODE_DIR[mode] / f"{tool}_time.txt"
    if not f.exists():
        return None
    vals = {}
    for line in open(f):
        k, _, v = line.strip().partition(",")
        vals[k] = v
    if "TotalTime" not in vals:
        return None
    return float(vals["TotalTime"]), "TotalTimeEstimated" in vals


def main():
    ap = argparse.ArgumentParser(description="Time against dataset size")
    ap.add_argument("--mode", default="perEdge", choices=["perEdge", "perBatch"])
    ap.add_argument("--datasets", nargs="+", required=True)
    ap.add_argument("--algorithms", nargs="+")
    ap.add_argument("--tools", nargs="*", default=["bubblegun", "bubblefinder"],
                    help="external tools: their per-edge time where it exists (hollow if estimated)")
    ap.add_argument("--estimate", nargs="*", default=[],
                    help="estimated times of runs that were not started, as ALGO:DATASET:HOURS (hollow)")
    ap.add_argument("--limit-days", type=float, default=2)
    ap.add_argument("--split", action="store_true",
                    help="two panels: our algorithms on a linear axis, the tools (and the best algorithm as "
                         "the reference) on a logarithmic one")
    ap.add_argument("--broken", action="store_true",
                    help="one plot with broken axes: small | medium and large, tools above | algorithms below")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    if not args.algorithms:
        args.algorithms = (["static", "incRecompute", "incNormalBST", "incBST", "incMinIBST", "incMinIBSTH1"]
                           if args.mode == "perEdge" else
                           ["static", "incRecomputeBatch", "incBSTBatch", "incBSTBatchH1", "incBSTBatchH2"])
    unit = 3600.0 if args.mode == "perEdge" else 1.0            # perEdge in hours, perBatch in seconds
    unit_name = "h" if args.mode == "perEdge" else "s"

    sizes, points = {}, {a: [] for a in args.algorithms}       # points: (dataset, edges, time, estimated)
    for d in args.datasets:
        folder = REPO / "output" / d / "timeFile" / MODE_DIR[args.mode]
        for a in args.algorithms:
            r = run_time(folder / f"{a}_time.txt")
            if r:
                sizes[d] = r[0] if args.mode == "perEdge" else (history_edges(d) or edgelist_edges(d) or r[0])
                points[a].append((d, sizes[d], r[1] / unit, False))
    if not sizes:
        sys.exit("no finished runs found")
    for t in args.tools:
        points[t] = []
        for d in args.datasets:
            r = tool_time(d, t, args.mode)
            if r and d in sizes:
                points[t].append((d, sizes[d], r[0] / unit, r[1]))
    for e in args.estimate:
        a, d, h = e.split(":")
        if d in sizes:
            points.setdefault(a, []).append((d, sizes[d], float(h) * 3600.0 / unit, True))

    plt.rcParams.update({"font.size": 13, "legend.fontsize": 11})
    limit = args.limit_days * 86400.0 / unit
    best = "incMinIBSTH1" if args.mode == "perEdge" else "incBSTBatchH2"

    def draw(ax, names, log):
        """the series in names on one axis; log: logarithmic y axis with the limit line"""
        for a in names:
            pts_all = points.get(a, [])
            if not pts_all:
                continue
            m, c = MARKERS.get(a, ("o", "black"))
            ms = 13 if m == "*" else 8
            name = TOOL_NAMES.get(a, SHORT_NAMES.get(a, a))
            first = True
            # measured runs of release 1 joined by a line; release 2 (another graph) and estimates single
            groups = [([p for p in pts_all if p[0].startswith("HPRC/") and not p[3]], True),
                      ([p for p in pts_all if not (p[0].startswith("HPRC/") and not p[3])], False)]
            for pts, joined in groups:
                pts = sorted(pts, key=lambda p: p[1])
                for p in (pts if not joined else [None]):
                    sel = pts if joined else [p]
                    if not sel:
                        continue
                    est = (not joined) and p[3]
                    ax.plot([q[1] for q in sel], [q[2] for q in sel], marker=m, color=c, markersize=ms,
                            markerfacecolor="white" if est else c, markeredgecolor=c, markeredgewidth=1.5,
                            linestyle=("-" if a == best else ":") if joined else "none", linewidth=1.5,
                            label=name if first else None)
                    first = False
        vals = [p[2] for a in names for p in points.get(a, [])]
        # dataset names next to the points of the best algorithm (it finishes everywhere), else of any series
        named = {}
        for a in [best] + list(names):
            for d, n, h, _ in points.get(a, []):
                if a in names or a == best:
                    named.setdefault(d, (n, h))
        for d, (n, h) in named.items():
            ax.annotate(label(d), (n, h), xytext=(0, -16), textcoords="offset points", ha="center", fontsize=10)
        ax.set_xscale("log")
        xs = [n for n in sizes.values()]
        ax.set_xlim(min(xs) * 0.92, max(xs) * 1.1)               # room for the dataset names at the edges
        ax.set_xlabel("inserted edges" + (" (updates)" if args.mode == "perEdge" else " (over all haplotypes)"))
        ax.grid(True, which="both", alpha=0.3, linestyle="--")
        if log:
            ax.axhline(limit, color="black", linestyle="--", linewidth=2)
            ax.text(max(sizes.values()), limit * 1.12, f"time limit ({args.limit_days:g} days)", fontsize=14,
                    fontweight="bold", ha="right")
            ax.set_yscale("log")
            lo, hi = min(vals) / 2.5, max(max(vals), limit) * 2
            ax.set_ylim(lo, hi)
            ticks = [mm * 10 ** e for e in range(-2, 7) for mm in (1, 2, 5) if lo <= mm * 10 ** e <= hi]
            ax.yaxis.set_major_locator(mticker.FixedLocator(ticks))
            ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f}" if v >= 1 else f"{v:g}"))
            ax.yaxis.set_minor_formatter(mticker.NullFormatter())
        else:
            ax.set_ylim(-0.1 * max(vals), max(vals) * 1.1)
            ax.set_yticks([t for t in ax.get_yticks() if t >= 0])

    def draw_broken():
        """one column per class of datasets, tools above and our algorithms below (y axis), own log range each.
        Only measured runs: an estimate counts as did not finish."""
        is_small = lambda d: d.startswith("Human/") or d == "Ecoli"
        is_large = lambda d: d.startswith("HPRC2/")
        if args.mode == "perBatch":      # medium on its own: the range of release 2 would flatten it
            cols = [("Small", is_small, 1.0, True), ("Medium (HPRC release 1)", lambda d: d.startswith("HPRC/"), 1.0, True),
                    ("Large (release 2)", is_large, 1.0, True)]
            widths = [1, 1.35, 0.42]
        else:                            # the tools did not finish on medium and large: no upper band there
            cols = [("Small", is_small, 1.0, True),
                    ("Medium (HPRC release 1) and large (release 2, chr22 v2)", lambda d: not is_small(d), unit, False)]
            widths = [1, 1.7]
        fig = plt.figure(figsize=(13.5 if args.mode == "perBatch" else 12, 8.8))
        gs = fig.add_gridspec(2, len(cols), height_ratios=[1, 2.2], width_ratios=widths, hspace=0.05, wspace=0.22)
        order = list(args.tools) + list(args.algorithms)
        drawn = set()
        for col, (title, in_col, cunit, has_upper) in enumerate(cols):
            if has_upper:
                bands = [(fig.add_subplot(gs[0, col]), True), (fig.add_subplot(gs[1, col]), False)]
            else:
                bands = [(fig.add_subplot(gs[:, col]), False)]
            for ax, upper_band in bands:
                vals = []
                for a in order:
                    m, c = MARKERS.get(a, ("o", "black"))
                    pts = [p for p in points.get(a, []) if in_col(p[0]) and not p[3] and (a in args.tools) == upper_band]
                    if not pts:
                        continue
                    drawn.add(a)
                    y = lambda p: p[2] * unit / cunit
                    vals += [y(p) for p in pts]
                    lines = {}
                    for p in sorted(pts, key=lambda p: p[1]):
                        lines.setdefault(p[0].split("/")[0], []).append(p)
                    for sel in lines.values():
                        # the best algorithm below the others: a smaller marker on top of it stays visible
                        ax.plot([p[1] for p in sel], [y(p) for p in sel], marker=m, color=c,
                                markersize=14 if m == "*" else 7, markeredgecolor=c, markeredgewidth=1.0, linewidth=1.5,
                                linestyle="-" if a == best else ":", zorder=2 if a == best else 3)
                if not upper_band and args.mode == "perEdge" and not has_upper:
                    lim = limit * unit / cunit
                    vals.append(lim)
                    ax.axhline(lim, color="black", linestyle="--", linewidth=2)
                    ax.text(max(sizes.values()) * 1.04, lim * 1.05, f"time limit ({args.limit_days:g} days)",
                            fontsize=12, fontweight="bold", ha="right")
                    ax.text(0.97, 0.04, "BubbleGun, BubbleFinder: did not finish on any dataset here\nOther missing points: did not finish",
                            transform=ax.transAxes, ha="right", va="bottom", fontsize=11,
                            bbox=dict(facecolor="white", edgecolor="black", boxstyle="round,pad=0.4"))
                ax.set_yscale("log")
                if vals:
                    wide = max(vals) / min(vals) >= 10
                    lo, hi = min(vals) / (1.5 if wide else 1.2), max(vals) * (1.5 if wide else 1.2)
                    ax.set_ylim(lo, hi)
                    steps = (1, 2, 5) if wide else (1, 1.5, 2, 3, 4, 5, 6, 8)
                    ticks = [mm * 10 ** e for e in range(-3, 8) for mm in steps if lo <= mm * 10 ** e <= hi]
                    ax.yaxis.set_major_locator(mticker.FixedLocator(ticks))
                ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f}" if v >= 1 else f"{v:g}"))
                ax.yaxis.set_minor_formatter(mticker.NullFormatter())
                ax.grid(True, which="major", alpha=0.35, linestyle="--")
                if not upper_band:
                    for d, n, h, _ in points.get(best, []):
                        if in_col(d):
                            ax.annotate(label(d), (n, h * unit / cunit), xytext=(0, -17),
                                        textcoords="offset points", ha="center", fontsize=9)
            bottom = bands[-1][0]
            xs = [n for d, n in sizes.items() if in_col(d)]
            for ax, _ in bands:
                ax.set_xscale("log")
                ax.set_xlim(min(xs) * 0.9, max(xs) * 1.1)
                ax.xaxis.set_minor_formatter(mticker.NullFormatter())
                ax.xaxis.set_minor_locator(mticker.NullLocator())
            if has_upper:                    # the break between the bands
                top = bands[0][0]
                top.spines["bottom"].set_visible(False)
                bottom.spines["top"].set_visible(False)
                top.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
                mark = dict(marker=[(-1, -0.6), (1, 0.6)], markersize=9, linestyle="none", color="black",
                            clip_on=False, mew=1.1)
                top.plot([0, 1], [0, 0], transform=top.transAxes, **mark)
                bottom.plot([0, 1], [1, 1], transform=bottom.transAxes, **mark)
            cand = ([7000, 10000, 20000] if is_small(next(d for d in sizes if in_col(d))) else
                    (xs if len(xs) == 1 else [1000000, 2000000, 3000000]))
            bottom.xaxis.set_major_locator(mticker.FixedLocator([x for x in cand if min(xs) * 0.9 <= x <= max(xs) * 1.1]))
            bottom.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
            bottom.set_xlabel("inserted edges" if len(xs) == 1 else
                              "inserted edges" + (" (updates)" if args.mode == "perEdge" else " (over all haplotypes)"))
            bottom.set_ylabel(f"total time ({'h' if cunit == 3600.0 else 's'}, log scale)")
            bands[0][0].set_title(title, fontsize=13)
        handles = []
        for a in order:
            if a in drawn:
                m, c = MARKERS.get(a, ("o", "black"))
                handles.append(plt.Line2D([], [], marker=m, color=c, markersize=14 if m == "*" else 7,
                                          linestyle="-" if a == best else ":",
                                          label="static" if a == "static" else TOOL_NAMES.get(a, SHORT_NAMES.get(a, a))))
        fig.legend(handles=handles, loc="upper center", ncol=len(handles), frameon=False, fontsize=11)
        fig.subplots_adjust(left=0.08, right=0.98, top=0.9, bottom=0.08)
        return fig

    subtitle = ("(filled: measured; hollow: estimated; no point: did not finish)"
                if any(p[3] for v in points.values() for p in v) else "(all points measured)")
    kind = "Per-edge" if args.mode == "perEdge" else "Per-batch"
    if args.broken:
        fig = draw_broken()
    elif args.split:
        # the algorithms on a linear axis (their differences are visible), the tools on a logarithmic one
        # (with the best algorithm as the reference)
        fig, (axl, axr) = plt.subplots(1, 2, figsize=(15, 6))
        draw(axl, args.algorithms, False)
        axl.set_ylabel(f"total time ({unit_name})")
        axl.set_title("Our algorithms (linear scale)", fontsize=13)
        axl.legend(loc="upper left")
        draw(axr, list(args.tools) + [best], True)
        axr.set_title("External tools (log scale)", fontsize=13)
        axr.legend(loc="upper left", bbox_to_anchor=(0.0, 0.86))
        fig.suptitle(f"{kind} algorithms on the HPRC chromosomes: time against size {subtitle}", fontsize=14)
        fig.tight_layout(rect=(0, 0, 1, 0.95))
    else:
        fig, ax = plt.subplots(figsize=(8.5, 6))
        draw(ax, list(points), True)
        ax.set_ylabel(f"total time ({unit_name})")
        ax.set_title(f"{kind} algorithms on the HPRC chromosomes: time against size\n{subtitle}", fontsize=13)
        ax.legend(loc="lower right" if args.mode == "perEdge" else "upper left", ncol=2)
        fig.tight_layout()
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    fig.savefig(args.out, dpi=200)
    print(f"Saved plot: {args.out}")
    for a, v in points.items():
        print(f"  {a}: " + ", ".join(f"{label(d)} {h:.1f} {unit_name}" + (" (est.)" if e else "")
                                     for d, _, h, e in sorted(v, key=lambda p: p[1])))


if __name__ == "__main__":
    main()
