r"""
plotcsv.py

Makes the performance plots of one dataset and mode from the CSV written
by parse_to_csv.py: time per update, time components, totals per
algorithm, and affected/explored vertices and edges.

=====================================================================
HOW TO RUN
=====================================================================

Run from the repository root, with the virtual environment active
(source venv/bin/activate), after Preprocessing/parse_to_csv.py.

    # E. coli, per-edge results, logarithmic y-axis
    python3 Preprocessing/plotcsv.py \
        --csv output/Ecoli/timeFile/perEdge/haplotype_build_comparison.csv \
        --out output/plots/perEdge/Ecoli --mode peredge --dataset Ecoli --scale log

    # chromosome 11, per-haplotype results, linear y-axis
    python3 Preprocessing/plotcsv.py \
        --csv output/Human/chromosome_11/timeFile/perbatch/haplotype_build_comparison.csv \
        --out output/plots/perBatch/chromosome_11 --mode perbatch --dataset chromosome_11

Options:

    --csv       haplotype_build_comparison.csv of one dataset and mode   (required)
    --out       folder for the PNG files (created if missing); by
                convention output/plots/<perEdge|perBatch>/<dataset>     (required)
    --mode      peredge | perbatch, must match the CSV                    (required)
    --dataset   name used in titles and file names, e.g. Ecoli, chr11     (required)
    --scale     linear | log y-axis                                       (default: linear)

=====================================================================
OUTPUT (PNG files in --out, named <dataset>_<mode>_<scale>_<plot>.png)
=====================================================================

    total_time_line       total time per update, one line per algorithm
    total_time_bar        total time over all updates, per algorithm
    components            time of each step (insert, DFS, topological
                          order, minI, maxO, superbubbles) per algorithm
    stacked_breakdown     the same steps stacked, for a few updates
    aggregated_versions   time of each step summed over all updates
    activity              part of the graph touched (affected / explored
                          vertices and edges): perbatch cumulative over
                          the genomes, edge-by-edge vs batch; peredge
                          distribution over the updates that reorder

Algorithms are shown with their display names (DISPLAY_NAMES below, e.g.
Inc-MinI-BST); the plotted algorithms are listed in ALGORITHMS below.
"""

import pandas as pd
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})
import sys
import os
import numpy as np # Import numpy for bar chart positioning
from matplotlib.lines import Line2D
import argparse

# Define the components and their corresponding column names

# all algorithms that can be plotted, in the order in which they were
# developed (legends, bars); each plot shows those present in the CSV
# (select_algorithms), e.g. incBSTBatch only in perBatch
ALL_ALGORITHMS = ["static", "incRecompute", "incNormalBST", "incBST", "incMinIBST", "incMinIBSTH1", "incBSTBatch"]
ALGORITHMS = list(ALL_ALGORITHMS)

# Names shown in legends (IDs and meanings: see ALGORITHMS.md)
DISPLAY_NAMES = {
    "static": "Static",
    "incRecompute": "Inc-Recompute",
    "incNormalBST": "Inc-Normal-BST",
    "incBST": "Inc-BST",
    "incMinIBST": "Inc-MinI-BST",
    "incMinIBSTH1": "Inc-MinI-BST-H1",
    "incBSTBatch": "Inc-BST-Batch",
}

# short labels under the bars of the stacked breakdown plots
SHORT_NAMES = {
    "static": "S",
    "incRecompute": "Rec",
    "incNormalBST": "NBST",
    "incBST": "BST",
    "incMinIBST": "MinI",
    "incMinIBSTH1": "H1",
    "incBSTBatch": "Batch",
}

COMPONENT_SUFFIX = {
    'UpdateDS': 'UpdateDSTime', 'DFS': 'DFSTime', 'TopologicalOrder': 'TopologicalOrderTime',
    'MinI': 'MinITime', 'MaxO': 'MaxOTime', 'SBCompute': 'SBComputeTime', 'Total': 'TotalTime',
}
# component -> one column per algorithm in ALGORITHMS (same order)
TIME_COMPONENTS = {comp: [a + suf for a in ALGORITHMS] for comp, suf in COMPONENT_SUFFIX.items()}


def select_algorithms(df):
    """Restrict ALGORITHMS / TIME_COMPONENTS to the algorithms in this CSV."""
    global ALGORITHMS, TIME_COMPONENTS
    ALGORITHMS = [a for a in ALL_ALGORITHMS if f"{a}TotalTime" in df.columns]
    TIME_COMPONENTS = {comp: [a + suf for a in ALGORITHMS] for comp, suf in COMPONENT_SUFFIX.items()}

# Define color palette for the stacked bar chart components
# --- MODIFIED: 'UpdateDS' color has been removed entirely ---
STACK_COLORS = {
    'UpdateDS': '#4C72B0',
    'DFS': '#55A868',
    'TopologicalOrder': '#DD8452',  
    'MinI': '#C44E52',
    'MaxO': '#8172B2',
    'SBCompute': '#CCB974',
}

# Define line styles for the three versions in line plots
VERSION_STYLES = {
    'static':  {'color': 'red',    'marker': 'o', 'linestyle': '-',   'linewidth': 0.7, 'markersize': 2},
    'incRecompute':     {'color': 'blue',   'marker': 's', 'linestyle': '--',  'linewidth': 0.7, 'markersize': 2},
    'incNormalBST':{'color': 'teal',  'marker': 'P', 'linestyle': (0, (3, 1, 1, 1)), 'linewidth': 0.7, 'markersize': 2},
    'incBST':     {'color': 'green',  'marker': '^', 'linestyle': '-.',  'linewidth': 0.7, 'markersize': 2},
    'incMinIBST':{'color': 'orange', 'marker': 'D', 'linestyle': ':',   'linewidth': 0.7, 'markersize': 2},
    'incMinIBSTH1':{'color': 'purple', 'marker': 'v', 'linestyle': '-',   'linewidth': 0.7, 'markersize': 2},
    'incBSTBatch':{'color': 'saddlebrown', 'marker': 'x', 'linestyle': '--',  'linewidth': 0.7, 'markersize': 2},
}

def build_plot_filename(dataset, mode, scale, plot_type):
    """
    Standardized filename generator.

    Example:
    GRCh38_chr1_peredge_log_components.png
    ecoli_perbatch_linear_total_bar.png
    """
    return f"{dataset}_{mode}_{scale}_{plot_type}.png"

# the edge-by-edge algorithm of the activity plot (all of them give the same
# affected / explored counts since they share the case logic)
EDGE_BY_EDGE_FOR_ACTIVITY = ["incMinIBSTH1", "incMinIBST", "incBST"]


def plot_activity(df, output_dir, dataset_name, mode, scale):
    """
    How much of the graph the algorithms touch: affected = reordered
    positions, explored = vertices reached by the searches, as a share of
    the graph (vertices: AffectedVertices / TotalVertices; edges, their
    out-degrees: TotalAffectedEdges / TotalEdges). Two panels, vertices and
    edges. One file: <dataset>_<mode>_<scale>_activity.png
      perbatch  cumulative over the genomes (absolute, log scale),
                edge-by-edge (the work of all its edges) and batch; the
                totals after the last genome are in the legend
      peredge   distribution over the edge updates that reorder (or search):
                most updates reorder nothing, the others often a large part
                of the graph; the counts are in the title
    """
    df_plot = df[df["Step"] != 0].copy()
    ebe = next((a for a in EDGE_BY_EDGE_FOR_ACTIVITY if f"{a}AffectedVertices" in df_plot.columns), None)
    algos = [a for a in [ebe] if a]
    # a count column with negative values overflowed in the program (32-bit
    # sum over a whole HPRC genome): such an algorithm is not plotted
    def overflowed(a):
        cols = [c for c in df_plot.columns if c.startswith(a) and c[len(a):] in
                ("AffectedVertices", "ExploredVertices", "TotalAffectedEdges", "TotalExploredEdges")]
        return any((pd.to_numeric(df_plot[c], errors="coerce") < 0).any() for c in cols)
    if ebe and overflowed(ebe):
        print(f"  [!] {dataset_name}: counts of {ebe} overflowed (negative), not plotted")
        algos = []
    if mode == "perbatch" and "incBSTBatchAffectedVertices" in df_plot.columns:
        algos.append("incBSTBatch")
    if not algos:
        return
    panels = [("Vertices", "AffectedVertices", "ExploredVertices", "TotalVertices"),
              ("Edges", "TotalAffectedEdges", "TotalExploredEdges", "TotalEdges")]

    def share(a, col, total):
        num = pd.to_numeric(df_plot.get(f"{a}{col}"), errors="coerce")
        den = pd.to_numeric(df_plot.get(f"{a}{total}"), errors="coerce")
        return 100.0 * num / den.where(den > 0)

    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    if mode == "perbatch":
        # cumulative: positions reordered (and vertices searched) up to genome
        # k; monotone, so genomes without any reorder leave no gap
        for ax, (what, aff, exp, tot) in zip(axes, panels):
            for a in algos:
                name = "Batch" if a == "incBSTBatch" else f"Edge-by-edge ({DISPLAY_NAMES.get(a, a)})"
                color = VERSION_STYLES[a]["color"]
                for col, kind, style in [(aff, "affected", "-"), (exp, "explored", ":")]:
                    vals = pd.to_numeric(df_plot.get(f"{a}{col}"), errors="coerce").fillna(0).cumsum()
                    ax.plot(df_plot["Step"], vals.where(vals > 0), color=color, linestyle=style, linewidth=1.4,
                            label=f"{kind}, {name}: {int(vals.iloc[-1]):,}" if len(vals) else f"{kind}, {name}")
            ax.set_yscale("log")
            ax.set_xlabel("Genome (step)")
            ax.set_ylabel(f"{what} touched, cumulative (log)")
            ax.set_title(what)
            ax.grid(True, which="both", linestyle="--", alpha=0.4)
            ax.legend(fontsize=9, title="total after the last genome", title_fontsize=9)
        fig.suptitle(f"Part of the graph touched, cumulative over the genomes — {dataset_name}\n"
                     "(affected: reordered positions; explored: vertices reached by the searches; "
                     "edges: their out-degrees)", fontsize=13)
    else:
        a = algos[0]
        n_updates = len(df_plot)
        n_touch = 0                      # updates that reorder (the same for every panel)

        def pct(x):
            return f"{x:.2g}%" if x < 10 else f"{x:.0f}%"
        for ax, (what, aff, exp, tot) in zip(axes, panels):
            aff_pct = share(a, aff, tot)
            exp_pct = share(a, exp, tot)
            aff_pct = aff_pct[aff_pct > 0].dropna()
            exp_pct = exp_pct[exp_pct > 0].dropna()
            values = pd.concat([aff_pct, exp_pct])
            if values.empty:
                ax.set_title(f"{what}: no update reorders")
                continue
            lo = max(values.min(), 1e-4)
            hi = max(values.max(), lo * 10)
            bins = np.logspace(np.log10(lo), np.log10(hi), 40)
            n_touch = max(n_touch, len(aff_pct), len(exp_pct))
            # legend: how much one such update touches (median, mean)
            ax.hist(aff_pct, bins=bins, alpha=0.6, color=VERSION_STYLES[a]["color"],
                    label=f"affected (positions repaired): median {pct(aff_pct.median())}, "
                          f"mean {pct(aff_pct.mean())}")
            ax.hist(exp_pct, bins=bins, alpha=0.6, color="gray",
                    label=f"explored (reached by the search): median {pct(exp_pct.median())}, "
                          f"mean {pct(exp_pct.mean())}")
            ax.set_xscale("log")
            ax.set_xlabel(f"% of the graph's {what.lower()} touched by one update (log)")
            ax.set_ylabel("Number of updates")
            ax.set_title(what)
            ax.legend(fontsize=9)
            ax.grid(True, which="both", linestyle="--", alpha=0.4)
        fig.suptitle(f"Size of the reorders of single edge updates — {dataset_name}, "
                     f"{DISPLAY_NAMES.get(a, a)}\n({n_touch:,} of {n_updates:,} updates reorder and are "
                     f"counted here; all other updates touch nothing)", fontsize=13)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, build_plot_filename(dataset_name, mode, scale, "activity")), dpi=200)
    plt.close(fig)


def plot_total_time_bar(df, output_dir, dataset_name, mode, scale):
    """
    Plot total time (sum across all steps excluding Base)
    for static, incRecompute, incBST as a single bar chart.
    """

    print("Generating Plot: Total Time Comparison Bar Chart...")

    df_plot = df[df["Step"] != 0].copy()

    versions = ALGORITHMS
    totals = []

    for v in versions:
        col = f"{v}TotalTime"
        if col not in df_plot.columns:
            raise KeyError(f"Missing column: {col}")
        totals.append(df_plot[col].sum())

    fig, ax = plt.subplots(figsize=(7, 5))

    # one colour per algorithm, as in the line plots; display names, slanted
    labels = [DISPLAY_NAMES.get(v, v) for v in versions]
    bars = ax.bar(labels, totals, color=[VERSION_STYLES[v]['color'] for v in versions], edgecolor="black", linewidth=0.6)
    plt.setp(ax.get_xticklabels(), rotation=20, ha="right")

    ax.set_ylabel("Total Time (ms)")
    if small_algorithms(dict(zip(versions, totals))):
        ax.set_yscale("log")
        ax.set_ylabel("Total Time (ms, log scale)")
    ax.set_title(f"Total Time Comparison — {dataset_name}")
    ax.grid(axis="y", linestyle="--", alpha=0.6)

    # Add value labels above bars
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2,
                height * 1.01,
                f"{height:.2f}",
                ha='center', va='bottom', fontsize=9)

    plt.tight_layout()

    filename = os.path.join(
    output_dir,
    build_plot_filename(dataset_name, mode, scale, "total_time_bar")
    )
    plt.savefig(filename, dpi=200, bbox_inches="tight")
    print(f"Saved total time bar plot to {filename}")
    plt.close(fig)


def clean_and_normalize_df(df):
    """
    1) Ensure Step is numeric and 'Base'->0
    2) Normalize expected column names (static/incRecompute/incBST + component)
    3) Coerce time columns to numeric and fill NaN with 0
    4) Fail fast with a readable error if a required column is missing
    """
    # 1) normalize Step
    if 'Step' in df.columns:
        df['Step'] = df['Step'].replace({'Base': 0})
        df['Step'] = pd.to_numeric(df['Step'], errors='coerce').fillna(0).astype(int)
    else:
        raise KeyError("Required column 'Step' not found in CSV")

    # 2) expected components and versions
    versions = ALGORITHMS
    # NOTE: 'UpdateDS' is kept here for data cleaning, but excluded from TIME_COMPONENTS for plotting
    components = ['UpdateDS', 'DFS', 'TopologicalOrder', 'MinI', 'MaxO', 'SBCompute', 'Total']  

    # 3) build canonical column names and check presence (allow small variations)
    # For robustness, accept lowercase variants or missing 'Time' suffix
    def find_col(df_cols, candidates):
        for c in candidates:
            if c in df_cols:
                return c
        return None

    cols = df.columns.tolist()
    for v in versions:
        for comp in components:
            # candidates: e.g., "staticUpdateDSTime", "staticUpdateDS", "simple_insert_time"
            candidates = [
                f"{v}{comp}Time",
                f"{v}{comp}",
                f"{v}_{comp}Time",
                f"{v}_{comp}",
                f"{v}{comp}time",
                f"{v}_{comp}_time"
            ]
            found = find_col(cols, candidates)
            if not found:
                # Not fatal for some components (you may allow missing SBCompute), but warn/raise as needed
                # Here we'll raise for core components, warn for SuperBubble
                if comp in ['UpdateDS', 'TopologicalOrder', 'DFS', 'MinI', 'MaxO']:
                    raise KeyError(f"Missing required column for plotting: tried {candidates}")
                else:
                    # Fill missing SBCompute/Total with zeros
                    df[f"{v}{comp}Time"] = 0.0
                    continue
            # rename to canonical form
            canon = f"{v}{comp}Time"
            if found != canon:
                df[canon] = df[found]
    # 4) coerce all canonical time columns to numeric
    for v in versions:
        for comp in components:
            canon = f"{v}{comp}Time"
            if canon in df.columns:
                df[canon] = pd.to_numeric(df[canon], errors='coerce').fillna(0.0)
    return df


def convert_seconds_to_ms(df):
    """Convert all canonical time columns from seconds to milliseconds."""
    df = df.copy()
    versions = ALGORITHMS
    components = ['UpdateDS', 'DFS', 'TopologicalOrder', 'MinI', 'MaxO', 'SBCompute', 'Total']

    for v in versions:
        for comp in components:
            col = f"{v}{comp}Time"
            if col in df.columns:
                df[col] = df[col] * 1000.0
    return df


def smooth_df_rolling(df, window=15):
    """
    Smooth time columns using rolling mean.
    Keeps Step unchanged; only smooths time columns.
    """
    df = df.sort_values("Step").copy()

    versions = ALGORITHMS
    components = ['UpdateDS','DFS','TopologicalOrder', 'MinI', 'MaxO', 'SBCompute', 'Total']  # UpdateDS excluded in plots

    for v in versions:
        for comp in components:
            col = f"{v}{comp}Time"
            if col in df.columns:
                df[col] = df[col].rolling(window=window, center=True, min_periods=1).mean()

    return df

def plot_separate_time_components(df, output_dir, dataset_name, mode, scale):
    """
    Generates a 3x2 grid of subplots showing time component performance
    (log scale) vs. Haplotype Step for all three versions.
    Includes a clear legend for the line colors.
    
    --- MODIFIED: Layout changed from 2x3 to 3x2 to accommodate 5 plots ---
    """
    print("Generating Plot 1: Time Component Comparison (Line Plot)...")
    
    # Filter out the 'Base' step which often has zero times and affects log scale visualization
    df_plot = df[df['Step'] != 0].copy()
    
    # Change from (2, 3) to (3, 2) for 5 plots
    fig, axes = plt.subplots(3, 2, figsize=(16, 12))
    axes = axes.flatten() # Flattens the 2D array of axes into a 1D array

    for ax in axes:
        if ax.get_legend() is not None:
            ax.get_legend().remove()

    i = 0
    # Loop through the components (DFS, MinI, MaxO, SBCompute, Total)
    for component, cols in list(TIME_COMPONENTS.items())[:-1]:
        
        ax = axes[i]
        
        # Plot each version
        eps = 1e-9  # prevents log-scale problems on zeros

        for v in ALGORITHMS:
            ax.plot(df_plot['Step'], df_plot[cols[ALGORITHMS.index(v)]] + eps, **VERSION_STYLES[v], label=v)

        ax.set_title(f'{component} Time vs. Haplotype Step — {dataset_name}', fontsize=12)
        if scale == "log":
            ax.set_yscale("log", nonpositive="clip")
        else:
            ax.set_yscale("linear")
        ylabel = f"{component} Time (ms)"
        if scale == "log":
            ylabel += " [Log Scale]"

        ax.set_ylabel(ylabel)
        ax.grid(True, which='both', ls='--', linewidth=0.5, alpha=0.7)
        
        # --- MODIFIED: X-label logic for 3x2 grid ---
        # Only add X-label to bottom row (indices 4 and 5 in the flattened array)
        if i >= 4:
            ax.set_xlabel('Haplotype Step (Ordinal)')
        
        i += 1
        
    # --- ADDED: Turn off the empty 6th subplot ---
    for j in range(i, len(axes)):
        axes[j].axis("off")
        
    
    legend_handles = [Line2D([0], [0], **VERSION_STYLES[v]) for v in ALGORITHMS]
    legend_labels = [DISPLAY_NAMES[a] for a in ALGORITHMS]


    # Reserve space at the top and add a single global legend
    fig.subplots_adjust(top=0.92)   # make room for legend
    fig.legend(
        legend_handles,
        legend_labels,
        loc='upper center',
        bbox_to_anchor=(0.5, 0.98),
        ncol=3,
        title='Algorithm (color)'
    )

    # tighten layout (respect the reserved top area) and save including legend
    plt.tight_layout(rect=[0, 0, 1, 0.9])


    plot_filename = os.path.join(
    output_dir,
    build_plot_filename(dataset_name, mode, scale, "components")
    )
    plt.savefig(plot_filename, dpi=200, bbox_inches='tight')
    print(f"Saved Plot 1 to {plot_filename}")
    plt.close(fig)  # Close the figure

def plot_total_time_line(df, output_dir, dataset_name, mode, scale):
    """
    Plot Total Time vs Step, one line per algorithm in the CSV (ALGORITHMS,
    set by select_algorithms). This shows the overall runtime behaviour.
    """

    print("Generating Plot: Total Time vs Step (Line Plot)...")

    df_plot = df[df["Step"] != 0].copy()

    fig, ax = plt.subplots(figsize=(10,6))

    eps = 1e-9

    for v in ALGORITHMS:
        col = f"{v}TotalTime"
        if col not in df_plot.columns:
            continue
        ax.plot(df_plot["Step"], df_plot[col] + eps, **VERSION_STYLES[v], label=DISPLAY_NAMES[v])

    if scale == "log":
        ax.set_yscale("log", nonpositive="clip")

    ax.set_xlabel("Step")
    ax.set_ylabel("Total Time (ms)")
    ax.set_title(f"Total Time vs Step — {dataset_name}")

    ax.grid(True, linestyle="--", alpha=0.6)
    ax.legend()

    plot_filename = os.path.join(
        output_dir,
        build_plot_filename(dataset_name, mode, scale, "total_time_line")
    )

    plt.tight_layout()
    plt.savefig(plot_filename, dpi=200)
    plt.close(fig)

    print(f"Saved Total Time Line Plot to {plot_filename}")


# an algorithm whose total is below 1/ZOOM_RATIO of the largest one is
# invisible in a linear stacked bar (e.g. static and Batch next to the
# edge-by-edge algorithms on HPRC perBatch): such plots get a second panel
# with only these algorithms on their own axis
ZOOM_RATIO = 20


def small_algorithms(totals):
    """the algorithms (of a {algorithm: total} dict) below 1/ZOOM_RATIO of the largest total"""
    top = max(totals.values()) if totals else 0
    small = [a for a, t in totals.items() if top > 0 and t < top / ZOOM_RATIO]
    return small if 0 < len(small) < len(totals) else []


def draw_stacked_groups(ax, df_plot, group_labels, group_axis_label, versions=None):
    """
    One group of stacked bars per row of df_plot (a genome or an edge
    update), one bar per algorithm, stacked by time component. The groups
    are 0.8 wide whatever the number of algorithms, centred on their tick,
    with a dotted line halfway between groups; each bar has its short name
    (SHORT_NAMES) under it and the groups their label below those. Returns
    the legend handles of the components.
    """
    components = ['UpdateDS', 'DFS', 'TopologicalOrder', 'MinI', 'MaxO', 'SBCompute']
    versions = list(versions or ALGORITHMS)
    n = len(versions)
    width = 0.8 / n
    groups = np.arange(len(df_plot))
    offsets = (np.arange(n) - (n - 1) / 2) * width
    bottoms = {v: np.zeros(len(df_plot)) for v in versions}

    handles = {}
    for comp in components:
        for k, v in enumerate(versions):
            col = f"{v}{comp}Time"
            vals = df_plot[col].fillna(0).values if col in df_plot.columns else np.zeros(len(df_plot))
            bars = ax.bar(groups + offsets[k], vals, width * 0.9, bottom=bottoms[v],
                          color=STACK_COLORS[comp], edgecolor="black", linewidth=0.4)
            handles.setdefault(comp, bars)
            bottoms[v] += vals

    # separators halfway between groups
    for x in groups[:-1] + 0.5:
        ax.axvline(x=x, color="gray", linestyle=":", linewidth=1.0, alpha=0.6)
    ax.set_xlim(-0.5, len(df_plot) - 0.5)

    # short algorithm names under every bar ...
    ax.set_xticks(np.concatenate([groups + o for o in offsets]))
    ax.set_xticklabels([SHORT_NAMES.get(v, v) for v in versions for _ in groups],
                       rotation=90, fontsize=8)
    ax.tick_params(axis="x", length=2, pad=2)
    # ... and the group labels below them
    sec = ax.secondary_xaxis(-0.09)
    sec.set_xticks(groups)
    sec.set_xticklabels(group_labels, fontsize=11)
    sec.tick_params(axis="x", length=0)
    sec.spines["bottom"].set_visible(False)
    sec.set_xlabel(group_axis_label, fontsize=13)
    return handles


def stacked_legend_and_note(ax, fig, handles):
    # outside the axes, right of the plot: inside it covers the last group
    ax.legend(list(handles.values()), list(handles.keys()), loc="upper left",
              bbox_to_anchor=(1.005, 1.0), title="Time component", fancybox=True)
    note = ",  ".join(f"{SHORT_NAMES.get(v, v)} = {DISPLAY_NAMES.get(v, v)}" for v in ALGORITHMS)
    fig.text(0.5, 0.005, note, ha="center", va="bottom", fontsize=9)

def plot_stacked_bar_breakdown(df, output_dir, dataset_name, mode, scale, num_steps=5):
    """
    Stacked bars of the time components of every algorithm for the first
    num_steps genomes (the first genomes change the graph the most, so they
    are the most expensive ones, not typical ones).
    """
    print(f"Generating Plot 2: Total Time Stacked Breakdown (First {num_steps} Steps)...")

    df_plot = df[df['Step'] != 0].head(num_steps).copy()
    if df_plot.empty:
        print("Not enough non-base data to generate the stacked bar chart.")
        return

    totals = {v: sum(df_plot[f"{v}{c}Time"].fillna(0).sum() for c in
                     ['UpdateDS', 'DFS', 'TopologicalOrder', 'MinI', 'MaxO', 'SBCompute']
                     if f"{v}{c}Time" in df_plot.columns) for v in ALGORITHMS}
    small = small_algorithms(totals)
    labels = df_plot['Step'].astype(str).tolist()
    if small:
        fig, (ax, axz) = plt.subplots(1, 2, figsize=(20, 8), gridspec_kw={"width_ratios": [3, 2]})
    else:
        fig, ax = plt.subplots(figsize=(15, 8))
    handles = draw_stacked_groups(ax, df_plot, labels, "Genome (step)")
    ax.set_title(f'Time per component for the first {num_steps} genomes — {dataset_name}', fontsize=16)
    ax.set_ylabel('Total Time (ms) [Sum of Components]')
    ax.grid(axis='y', ls='--', alpha=0.6)
    if small:
        draw_stacked_groups(axz, df_plot, labels, "Genome (step)", versions=small)
        axz.set_title(f"Zoom: {', '.join(SHORT_NAMES.get(v, v) for v in small)} (own scale; below 1/{ZOOM_RATIO} of the largest)",
                      fontsize=12)
        axz.grid(axis='y', ls='--', alpha=0.6)
        stacked_legend_and_note(axz, fig, handles)
    else:
        stacked_legend_and_note(ax, fig, handles)
    plt.tight_layout(rect=(0, 0.03, 1, 1))

    plot_filename = os.path.join(
        output_dir,
        build_plot_filename(dataset_name, mode, scale, "stacked_breakdown")
    )
    plt.savefig(plot_filename, dpi=200, bbox_inches="tight")
    print(f"Saved Plot 2 to {plot_filename}")
    plt.close(fig)



def plot_aggregated_version_bars(df, output_dir, dataset_name, mode, scale):
    """
    Produce a single stacked bar per version (static / incRecompute / incBST).
    Each bar is the sum over all steps (excluding Step==0/'Base') of the component times.
    This yields exactly three bars and a clear component legend.
    """
    print("Generating Plot: Aggregated stacked bars (one bar per version)...")
    df_plot = df[df['Step'] != 0].copy()
    if df_plot.empty:
        print("No non-base rows found; skipping aggregated bar plot.")
        return

    versions = ALGORITHMS
    # --- MODIFIED: 'UpdateDS' component has been removed ---
    components = ['UpdateDS','DFS','TopologicalOrder', 'MinI', 'MaxO', 'SBCompute']

    # compute sums for each (version, component)
    sums = {v: [] for v in versions}
    for comp in components:
        for v in versions:
            col = f'{v}{comp}Time'
            if col not in df_plot.columns:
                raise KeyError(f"Expected column '{col}' not found.")
            sums[v].append(df_plot[col].sum())

    # plotting (with a zoom panel for algorithms that would be invisible)
    small = small_algorithms({v: sum(sums[v]) for v in versions})
    if small:
        fig, (ax, axz) = plt.subplots(1, 2, figsize=(13, 6), gridspec_kw={"width_ratios": [3, 2]})
    else:
        fig, ax = plt.subplots(figsize=(8, 6))

    def stacked(axis, vs):
        x = np.arange(len(vs))
        bottom = np.zeros(len(vs))
        for i, comp in enumerate(components):
            values = [sums[v][i] for v in vs]
            axis.bar(x, values, bottom=bottom, width=0.6, label=comp, color=STACK_COLORS[comp],
                     edgecolor='black', linewidth=0.4)
            bottom += values
        axis.set_xticks(x)
        axis.set_xticklabels([DISPLAY_NAMES.get(a, a) for a in vs], rotation=20, ha="right")
        axis.grid(axis='y', ls='--', alpha=0.5)

    stacked(ax, list(versions))
    ax.set_ylabel('Time (ms) — summed across steps')
    ax.set_title(f'Aggregated Total Time Breakdown by Version --{dataset_name}')
    if small:
        stacked(axz, small)
        axz.set_title(f"Zoom (own scale; below 1/{ZOOM_RATIO} of the largest)", fontsize=10)
        axz.legend(title='Component', loc='upper left', bbox_to_anchor=(1.01, 1.0))
    else:
        ax.legend(title='Component', loc='upper right')
    plt.tight_layout()


    plot_filename = os.path.join(
        output_dir,
        build_plot_filename(dataset_name, mode, scale, "aggregated_versions")
    )
    plt.savefig(plot_filename, dpi=200)
    print(f"Saved aggregated bar plot to {plot_filename}")
    plt.close(fig)


def prepare_peredge_df(df, max_points=5000):
    # Remove Base row
    df = df[df["Step"] != 0].copy()

    # Step already global integer, but enforce
    df["Step"] = pd.to_numeric(df["Step"], errors="coerce").fillna(0).astype(int)
    df = df[df["Step"] > 0].copy()

    # Downsample uniformly if too large (for speed + file size)
    if len(df) > max_points:
        idx = np.linspace(0, len(df) - 1, max_points).astype(int)
        df = df.iloc[idx].copy()

    return df


def plot_separate_time_components_peredge(df, output_dir, dataset_name, mode, scale):
    print("Generating Plot 1 (perEdge): Time Component Comparison (Line Plot)...")

    df_plot = df.copy()

    fig, axes = plt.subplots(3, 2, figsize=(16, 12))
    axes = axes.flatten()

    for ax in axes:
        if ax.get_legend() is not None:
            ax.get_legend().remove()

    i = 0
    for component, cols in list(TIME_COMPONENTS.items())[:-1]:
        if i >= 6:
            break

        ax = axes[i]
        eps = 1e-9

        for v in ALGORITHMS:
            ax.plot(df_plot['Step'], df_plot[cols[ALGORITHMS.index(v)]] + eps, **VERSION_STYLES[v], label=v)

        ax.set_title(f"{component} Time vs Edge Updates — {dataset_name}", fontsize=12)
        if scale == "log":
            ax.set_yscale("log", nonpositive="clip")
        else:
            ax.set_yscale("linear")
        ylabel = f"{component} Time (ms)"
        if scale == "log":
            ylabel += " [Log Scale]"  
        ax.set_ylabel(ylabel)
        ax.grid(True, which="both", ls="--", linewidth=0.5, alpha=0.7)

        if i >= 4:
            ax.set_xlabel("Edge Update ID (Global Step)")

        i += 1

    for j in range(i, len(axes)):
        axes[j].axis("off")

    handles = [Line2D([0], [0], **VERSION_STYLES[v]) for v in ALGORITHMS]
    labels = [DISPLAY_NAMES[a] for a in ALGORITHMS]

    fig.subplots_adjust(top=0.92)
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.98), ncol=3)

    plt.tight_layout(rect=[0, 0, 1, 0.9])
    


    plot_filename = os.path.join(
        output_dir,
        build_plot_filename(dataset_name, mode, scale, "components")
    )
    
    plt.savefig(plot_filename, dpi=200, bbox_inches="tight")
    print(f"Saved Plot 1 to {plot_filename}")
    plt.close(fig)


def plot_stacked_bar_breakdown_peredge(df, output_dir, dataset_name, mode, scale, num_steps=10):
    """Like plot_stacked_bar_breakdown, for the first num_steps edge updates."""
    print(f"Generating Plot 2 (perEdge): Total Time Stacked Breakdown (First {num_steps} Edge Updates)...")

    df_plot = df.head(num_steps).copy()
    if df_plot.empty:
        print("Not enough data to generate stacked bar chart.")
        return

    fig, ax = plt.subplots(figsize=(15, 8))
    handles = draw_stacked_groups(ax, df_plot, df_plot["Step"].astype(str).tolist(), "Edge update (global step)")
    ax.set_title(f"Time per component for the first {num_steps} edge updates — {dataset_name}", fontsize=16)
    ax.set_ylabel("Time (ms)")
    ax.grid(True, axis="y", linestyle="--", alpha=0.5)
    stacked_legend_and_note(ax, fig, handles)
    plt.tight_layout(rect=(0, 0.03, 1, 1))

    plot_filename = os.path.join(
        output_dir,
        build_plot_filename(dataset_name, mode, scale, "stacked_breakdown")
    )
    plt.savefig(plot_filename, dpi=200, bbox_inches="tight")
    print(f"Saved Plot 2 to {plot_filename}")
    plt.close(fig)

def run_perbatch(csv_path, output_dir,dataset_name, scale):
    df = pd.read_csv(csv_path)
    select_algorithms(df)          # before cleaning: it checks the columns of these algorithms only
    df = clean_and_normalize_df(df)
    # convert to ms
    df = convert_seconds_to_ms(df)

    # smoothing (batch steps are fewer, so small window)
    df_smooth = smooth_df_rolling(df, window=3)
    mode = "perbatch"

    plot_separate_time_components(df_smooth, output_dir, dataset_name, mode, scale)

    plot_stacked_bar_breakdown(df, output_dir, dataset_name, mode, scale, num_steps=5)

    plot_aggregated_version_bars(df, output_dir, dataset_name, mode, scale)

    plot_total_time_bar(df, output_dir, dataset_name, mode, scale)

    plot_total_time_line(df_smooth, output_dir, dataset_name, mode, scale)

    plot_activity(df, output_dir, dataset_name, mode, scale)

def plot_chromosome_comparison(csv_files, labels, output_dir):

    # the algorithms of the first CSV (all CSVs of one call have the same mode)
    select_algorithms(pd.read_csv(csv_files[0], nrows=1))
    totals = {alg: [] for alg in ALGORITHMS}

    for file in csv_files:

        df = pd.read_csv(file)
        df = clean_and_normalize_df(df)
        df = convert_seconds_to_ms(df)

        df = df[df["Step"] != 0]

        for alg in ALGORITHMS:
            col = f"{alg}TotalTime"
            totals[alg].append(df[col].sum())

    x = np.arange(len(labels))
    width = 0.8 / len(ALGORITHMS)

    fig, ax = plt.subplots(figsize=(10,6))

    for i, alg in enumerate(ALGORITHMS):

        ax.bar(
            x + i*width,
            totals[alg],
            width,
            label=DISPLAY_NAMES.get(alg, alg),
            color=VERSION_STYLES[alg]['color']
        )

    ax.set_xticks(x + width * (len(ALGORITHMS) - 1) / 2)
    ax.set_xticklabels(labels)

    ax.set_ylabel("Total Time (ms)")
    ax.set_title("Chromosome Algorithm Comparison")

    ax.legend()
    ax.grid(axis="y", linestyle="--", alpha=0.6)

    plt.tight_layout()

    filename = os.path.join(output_dir, "chromosome_comparison.png")
    plt.savefig(filename, dpi=200)
    plt.close(fig)

def run_peredge(csv_path, output_dir,dataset_name, scale):
    df = pd.read_csv(csv_path)
    select_algorithms(df)          # before cleaning: it checks the columns of these algorithms only
    df = clean_and_normalize_df(df)

    # remove Base and sort
    df = df[df["Step"] != 0].copy()
    df = df.sort_values("Step")

    # convert to ms
    df = convert_seconds_to_ms(df)

    # Thinning (downsampling) to at most max_points evenly spread rows, ONLY
    # for the line plots (a curve over a million updates is slow to draw and
    # looks the same with 5,000 points). The bar plots add the times up and
    # the stacked plot shows the first updates: they get every row. (They
    # used to get the thinned rows too: chr22 bars showed 9.48 s of static's
    # 14.45 s, and the stacked plot skipped every third update.)
    max_points = 5000
    def thin(d):
        if len(d) <= max_points:
            return d
        idx = np.linspace(0, len(d) - 1, max_points).astype(int)
        return d.iloc[idx].copy()

    # smoothing over all updates (bigger window), then thinning
    df_smooth = thin(smooth_df_rolling(df, window=50))

    mode = "peredge"

    # line plots: thinned
    plot_separate_time_components_peredge(df_smooth, output_dir, dataset_name, mode, scale)
    plot_total_time_line(df_smooth, output_dir, dataset_name, mode, scale)

    # bars, the first updates and the distribution of the reorders: every row
    plot_activity(df, output_dir, dataset_name, mode, scale)
    plot_stacked_bar_breakdown_peredge(df, output_dir, dataset_name, mode, scale, num_steps=10)

    plot_aggregated_version_bars(df, output_dir, dataset_name, mode, scale)

    plot_total_time_bar(df, output_dir, dataset_name, mode, scale)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True, help="Path to CSV file generated by parse_to_csv.py")
    parser.add_argument("--out", required=True, help="Output directory for plots")
    parser.add_argument("--mode", choices=["perbatch", "peredge"], required=True)
    parser.add_argument("--dataset", required=True, help="Dataset name")
    parser.add_argument(
    "--scale",
    choices=["linear", "log"],
    default="linear",
    help="Y-axis scale for line plots"
    )


    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    if args.mode == "perbatch":
        run_perbatch(args.csv, args.out, args.dataset, args.scale)
    else:
        run_peredge(args.csv, args.out,  args.dataset, args.scale)

    print("✅ Plot generation complete!")

if __name__ == "__main__":

    main()
    """
    Main entry point for the script.
    Expects two command-line arguments:
    1. The path to the input CSV file.
    2. The path to the output directory for plots.
    
    Example:

    python plotcsv.py --csv chr5.csv --out plots --mode peredge --dataset GRCh38_chr5 --scale log

    python plotcsv.py --csv chr5.csv --out plots --mode perbatch --dataset GRCh38_chr5 --scale log
    python plotcsv.py --csv chr5.csv --out plots --mode peredge --dataset GRCh38_chr5 --scale linear
    """
    
    # Ensure the user has the required packages
    # You can install them with:
    # pip install pandas matplotlib numpy
    """
    # Check if a filename AND output dir were provided as command-line arguments
    if len(sys.argv) > 2:
        # sys.argv[0] is the script name, sys.argv[1] is the input file, sys.argv[2] is the output dir
        filepath = sys.argv[1]
        output_path = sys.argv[2]
        plot_haplotype_data(filepath, output_path)
    else:
        # If not enough arguments are provided, show the user how to run the script
        print("Usage: python plotcsv.py <path_to_your_csv_file> <path_to_output_directory>")
        print(f"Example: python plotcsv.py haplotype_build_comparison.csv plots/")
    """    