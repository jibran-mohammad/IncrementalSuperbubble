#!/usr/bin/env python3
"""
groupedTable.py - the long time table of the paper (every step of every
algorithm per dataset, with the speedup over static) with the datasets in
three groups: small (human chr11, chr17, chr18, chr22, chrX), medium (the same
chromosomes of HPRC release 1) and large (HPRC release 2 chr22). Same layout
as the long table of the paper: one block per step (Insert = update of the
data structures, DFS, Top. Order, MinI, MaxO, SBCompute, Total); in a block
static first, then every other algorithm with its time and, below it, its
speedup over static. A Total column per group (small, medium): the sum over
its datasets (DNF if one of them did not finish); the large group has one
dataset, a Total would repeat it.

Also the tables of the explored / affected vertices and edges (one for the
vertices, one for the edges; perEdge and perBatch in blocks), compared with
what static processes, in the layout of the paper:
    <outdir>/table_vertex_processed.tex, <outdir>/table_edge_processed.tex

And the tables of the external tools (BubbleGun, BubbleFinder) against our
algorithms, one per mode, BubbleGun as the reference: BubbleGun, BubbleFinder,
then static and the other algorithms, each with its speedup over BubbleGun:
    <outdir>/table_tools_edge.tex, <outdir>/table_tools_batch.tex
Our algorithms: the wall time of the whole run (the BENCHMARK runs, measured
like the tools); where a dataset has none, the sum of the timed steps of the
normal run, marked with an asterisk.

Writes two tables, one per mode:
    <outdir>/table_full_time_edge.tex    perEdge: static, Rec, MinI, H1
    <outdir>/table_full_time_batch.tex   perBatch: static, H1, BRec, Batch, BatchH1
and one small table per step and mode (Insert, DFS, Top. Order, MinI, MaxO,
SBCompute, Total), the same layout and groups: static against the final algorithm
(perEdge: H1, perBatch: BatchH1), with its speedup:
    <outdir>/table_<step>_edge.tex       e.g. table_dfs_edge.tex
    <outdir>/table_<step>_batch.tex
The MaxO tables also show the data structure: perEdge static, BST and the
faster of MinI and H1 (maxO by comparison); perBatch static, Batch (BST) and
BatchH1.

Run from the repository root:
    python3 Preprocessing/groupedTable.py --outdir output/paperTables/grouped


Times are the per-update timers of the normal build summed over all updates
(as in paperTables.py); a run that did not finish, was not started or took
more than the time limit is DNF, with no speedup.
"""
import argparse
import datetime
import os
import pathlib
import sys

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))
from Preprocessing.csvToLatex import (ALGORITHMS, EDGE_BY_EDGE_FOR_ACTIVITY, SHORT_NAMES,   # noqa: E402
                                      TIME_NOTE, abbreviation_note, big_number, normalize_dataframe,
                                      ratio_text, sum_metric)
from Preprocessing.paperTables import (DNF, DNF_NOTE, LIMIT_S, STEPS, fmt, ratio,        # noqa: E402
                                       run_finished)
from Preprocessing.toolComparisonToLatex import (TOOL_LABELS, read_pipeline_time,      # noqa: E402
                                                 read_tool_time)
from pipelineRunner import TOOL_NAME as PIPELINE_DIR                                    # noqa: E402

CHROMS = ["11", "17", "18", "22", "X"]
GROUPS = [("Small", [(f"chr{c}", f"Human/chromosome_{c}") for c in CHROMS], True),
          ("Medium", [(f"chr{c}", f"HPRC/chromosome_{c}") for c in CHROMS], True),
          ("Large", [("chr22 v2", "HPRC2/chromosome_22")], False)]
MODES = {"perEdge": ("perEdge", ["static", "incRecompute", "incMinIBST", "incMinIBSTH1"], "edge"),
         "perBatch": ("perbatch", ["static", "incMinIBSTH1", "incRecomputeBatch", "incBSTBatch", "incBSTBatchH1", "incBSTBatchH2"],
                      "batch")}
STEP_LABEL = {"UpdateDS": "Insert", "DFS": "DFS", "TopoOrder": "Top. Order", "MinI": "MinI", "MaxO": "MaxO",
              "SBCompute": "SBCompute", "Total": "Total"}


def name_of(a):
    return "static" if a == "static" else SHORT_NAMES.get(a, a)


STEP_SLUG = {"UpdateDS": "insert", "DFS": "dfs", "TopoOrder": "toporder", "MinI": "mini", "MaxO": "maxo",
             "SBCompute": "sbcompute", "Total": "total"}
FINAL = {"perEdge": "incMinIBSTH1", "perBatch": "incBSTBatchH1"}
BST_VARIANT = {"perEdge": "incBST", "perBatch": "incBSTBatch"}


def load(mode):
    """the columns of a mode: (group, label, csv path, DataFrame or None)"""
    mode_dir = MODES[mode][0]
    cols = []
    for group, members, _ in GROUPS:
        for label, dataset in members:
            path = REPO / "output" / dataset / "timeFile" / mode_dir / "haplotype_build_comparison.csv"
            df = normalize_dataframe(pd.read_csv(path, low_memory=False), str(path)) if path.exists() else None
            cols.append((group, label, path, df))
            print(f"  {mode} {label} ({group})")
    return cols


def value(a, c, suffix):
    group, label, path, df = c
    if df is None or not run_finished(str(path), a, True):
        return DNF
    run = sum_metric(df, a, "TotalTime")
    if not pd.isna(run) and run > LIMIT_S:
        return DNF
    return sum_metric(df, a, suffix)


def group_total(vals):
    if DNF in vals:
        return DNF
    vals = [v for v in vals if not pd.isna(v)]
    return sum(vals) if vals else float("nan")


def layout_of(cols, totals=True):
    """the columns of the table: datasets and, after the datasets of a group, its Total (totals)"""
    layout = []                                    # (kind, group, column object or None)
    for group, members, with_total in GROUPS:
        for c in cols:
            if c[0] == group:
                layout.append(("data", group, c))
        if with_total and totals:
            layout.append(("total", group, None))
    return layout


def row_values(layout, a, suffix):
    out, per_group = [], {}
    for kind, group, c in layout:
        if kind == "data":
            v = value(a, c, suffix)
            per_group.setdefault(group, []).append(v)
            out.append(v)
        else:
            out.append(group_total(per_group[group]))
    return out


def head_lines(layout, box="adjustbox", totals=True):
    """\\begin{table} ... the group heading and the dataset heading rows; box = adjustbox | resizebox"""
    spec = "|l|"
    for group, members, with_total in GROUPS:
        k = len(members) + (1 if with_total and totals else 0)
        spec += "|".join(["c"] * k) + ("||" if group != GROUPS[-1][0] else "|")
    opening = (r"\begin{adjustbox}{max width=\linewidth}" if box == "adjustbox"
               else r"\resizebox{\textwidth}{!}{%")
    lines = [r"\begin{table}[ht]", r"\centering", opening, rf"\begin{{tabular}}{{{spec}}}", r"\hline"]
    heads = []
    for group, members, with_total in GROUPS:
        k = len(members) + (1 if with_total and totals else 0)
        end = "|" if group == GROUPS[-1][0] else "||"
        heads.append(rf"\multicolumn{{{k}}}{{c{end}}}{{\textbf{{{group}}}}}")
    lines += [" & " + " & ".join(heads) + r" \\", r"\hline"]
    names = [r"\textbf{" + (c[1] if kind == "data" else "Total") + "}" for kind, g, c in layout]
    lines += [r"\textbf{Dataset} & " + " & ".join(names) + r" \\", r"\hline\hline" if box == "adjustbox" else r"\hline"]
    return lines


def build_steps(mode, cols, outdir):
    """one small table per step: static against the final algorithm (MaxO: also the BST variant)"""
    mode_dir, algos, tag = MODES[mode]
    layout = layout_of(cols)
    final, bst = FINAL[mode], BST_VARIANT[mode]
    mode_text = "edge updates (perEdge)" if mode == "perEdge" else "haplotypes (perBatch)"
    for step, suffix in STEPS:
        base = row_values(layout, "static", suffix)
        if step == "MaxO":
            if mode == "perEdge":
                # MinI and H1 compute maxO by comparison: the faster of the two, over the
                # datasets where both finished
                mi, h1 = row_values(layout, "incMinIBST", suffix), row_values(layout, "incMinIBSTH1", suffix)
                both = [(x, y) for (k, _, _), x, y in zip(layout, mi, h1)
                        if k == "data" and DNF not in (x, y) and not pd.isna(x) and not pd.isna(y)]
                fast = "incMinIBST" if sum(x for x, _ in both) < sum(y for _, y in both) else "incMinIBSTH1"
                print(f"  perEdge MaxO: {SHORT_NAMES[fast]} is the faster of MinI / H1 "
                      f"({sum(x for x, _ in both):.4f} s vs {sum(y for _, y in both):.4f} s)")
                shown = [bst, fast]
            else:
                shown = [bst, final]
        else:
            shown = [final]
        lines = head_lines(layout)
        lines.append(r"\texttt{static} & " + " & ".join(fmt(v, 4) for v in base) + r" \\")
        lines.append(r"\hline")
        for a in shown:
            vals = row_values(layout, a, suffix)
            lines.append(rf"\texttt{{{name_of(a)}}} & " + " & ".join(fmt(v, 4) for v in vals) + r" \\")
            lines.append("Speedup & " + " & ".join(fmt(ratio(b, v), 2) for b, v in zip(base, vals)) + r" \\")
            lines.append(r"\hline")
        extra = ""
        if step == "MaxO":
            extra = (" BST keeps maxO in a tree, H1 computes it by comparison." if mode == "perEdge" else
                     " Batch keeps maxO in a tree.")
        caption = (f"{STEP_LABEL[step]} time (s) for {mode_text}, small / medium / large datasets; speedup = "
                   f"static / algorithm; Total = sum over the group; DNF = did not finish within 2 days.{extra}")
        lines += [r"\end{tabular}", r"\end{adjustbox}", rf"\caption{{{caption}}}",
                  rf"\label{{tab:{STEP_SLUG[step]}_{tag}}}", r"\end{table}"]
        out = os.path.join(outdir, f"table_{STEP_SLUG[step]}_{tag}.tex")
        with open(out, "w") as f:
            f.write("\n".join(lines) + "\n")
        print(f"wrote {out}")


ACT = {"V": ("ExploredVertices", "AffectedVertices", "TotalVertices", "vertices", "V", "vertex"),
       "E": ("TotalExploredEdges", "TotalAffectedEdges", "TotalEdges", "edges", "E", "edge")}


def read_activity(path, mode, what):
    """(explored, affected, static) summed over all updates of one CSV, None if not available.
    explored / affected: the edge-by-edge algorithm (perEdge: H1) or the batch algorithm (perBatch);
    static: the graph's size summed over the updates - a static run does a full DFS of the whole
    graph after every update, so this is what it processes (logged by the incremental runs).
    A negative count (a 32-bit overflow in the program) gives None."""
    if not path.exists():
        return None
    explored_col, affected_col, size_col = ACT[what][:3]
    header = list(pd.read_csv(path, nrows=0).columns)
    if mode == "perBatch":
        algo = "incBSTBatch" if f"incBSTBatch{explored_col}" in header else None
    else:
        algo = next((a for a in EDGE_BY_EDGE_FOR_ACTIVITY if f"{a}{explored_col}" in header), None)
    size_algo = next((a for a in ALGORITHMS if f"{a}{size_col}" in header), None)
    if algo is None or size_algo is None:
        return None
    use = (["Step"] if "Step" in header else []) + [f"{algo}{explored_col}", f"{algo}{affected_col}",
                                                    f"{size_algo}{size_col}"]
    df = pd.read_csv(path, usecols=use, low_memory=False)
    if "Step" in df.columns:
        df = df[df["Step"].astype(str) != "Base"]

    def total(col):
        v = pd.to_numeric(df[col], errors="coerce").fillna(0)
        return None if (v < 0).any() else float(v.sum())
    return total(f"{algo}{explored_col}"), total(f"{algo}{affected_col}"), total(f"{size_algo}{size_col}")


def build_activity(cols_by_mode, outdir):
    """the tables of the explored / affected vertices and edges: perEdge and perBatch in blocks"""
    for what in ("V", "E"):
        word, letter, fname = ACT[what][3], ACT[what][4], ACT[what][5]
        layout = layout_of(cols_by_mode["perEdge"])
        lines = head_lines(layout, "resizebox")
        n = len(layout)
        for bi, mode in enumerate(("perEdge", "perBatch")):
            data = {}                              # column index -> (explored, affected, static)
            for i, (kind, group, c) in enumerate(layout):
                if kind == "data":
                    path = next(x[2] for x in cols_by_mode[mode] if x[1] == c[1] and x[0] == group)
                    data[i] = read_activity(path, mode, what)
            per_group = {}
            for i, (kind, group, c) in enumerate(layout):
                if kind == "data":
                    per_group.setdefault(group, []).append(data[i])
                else:
                    ds = per_group[group]
                    data[i] = None if any(d is None or None in d for d in ds) else tuple(
                        sum(d[k] for d in ds) for k in range(3))

            def cell(i, k):
                d = data[i]
                return "-" if d is None or d[k] is None else big_number(d[k])

            def ratio_cell(i, a, b):
                d = data[i]
                return "-" if d is None or None in (d[a], d[b]) or d[b] == 0 else ratio_text(d[a] / d[b])
            order = range(n)
            lines.append(rf"\multicolumn{{{n + 1}}}{{|l|}}{{\textbf{{{mode}}}}} \\")
            lines.append(r"\hline")
            lines.append("Explored & " + " & ".join(cell(i, 0) for i in order) + r" \\")
            lines.append(r"\hline")
            lines.append("Affected & " + " & ".join(cell(i, 1) for i in order) + r" \\")
            lines.append("Ratio & " + " & ".join(ratio_cell(i, 1, 0) for i in order) + r" \\")
            lines.append(r"\hline")
            lines.append(rf"Static ($\sum$ {letter}) & " + " & ".join(cell(i, 2) for i in order) + r" \\")
            lines.append("Ratio & " + " & ".join(ratio_cell(i, 2, 0) for i in order) + r" \\")
            lines.append(r"\hline")
        caption = (f"Total {word} processed, small / medium / large datasets, perEdge and perBatch: Explored = "
                   f"reached by the searches, Affected = positions reordered; Static = the graph's {word} summed "
                   "over the updates, i.e. what static processes; Ratio = relative to Explored (Affected / "
                   "Explored; Static / Explored).")
        lines += [r"\end{tabular}%", r"}", rf"\caption{{{caption}}}", rf"\label{{tab:{fname}_processed}}",
                  r"\end{table}"]
        out = os.path.join(outdir, f"table_{fname}_processed.tex")
        with open(out, "w") as f:
            f.write("\n".join(lines) + "\n")
        print(f"wrote {out}")


def memory_value(c, a):
    """peak memory of an algorithm on a dataset (MB): the peak of the program minus its own memory before
    the input is read (peak_memory.csv next to the CSV); DNF if the run did not finish"""
    group, label, path, df = c
    if not run_finished(str(path), a, True):
        return DNF
    f = path.parent / "peak_memory.csv"
    if not f.exists():
        return float("nan")
    m = pd.read_csv(f)
    r = m[m["Algorithm"] == a]
    if r.empty:
        return float("nan")
    base = r["BaselineMemoryKB"].iloc[0] if "BaselineMemoryKB" in m.columns and pd.notna(r["BaselineMemoryKB"].iloc[0]) else 0
    return (r["PeakMemoryKB"].iloc[0] - base) / 1024.0


def build_memory(cols_by_mode, outdir):
    """one table per mode: the peak memory of every algorithm, with its ratio to static; no Total column"""
    for mode in ("perEdge", "perBatch"):
        cols = cols_by_mode[mode]
        mode_dir, algos, tag = MODES[mode]
        layout = layout_of(cols, totals=False)
        lines = head_lines(layout, totals=False)
        base = [memory_value(c, "static") for _, _, c in layout]
        for ai, a in enumerate(algos):
            vals = [memory_value(c, a) for _, _, c in layout]
            lines.append(rf"\texttt{{{name_of(a)}}} & " + " & ".join(fmt(v, 1) for v in vals) + r" \\")
            if a != "static":
                lines.append("Ratio & " + " & ".join(
                    fmt(None if DNF in (b, v) or pd.isna(b) or pd.isna(v) else v / b, 2)
                    for b, v in zip(base, vals)) + r" \\")
            lines.append(r"\hline")
        mode_text = "edge updates (perEdge)" if mode == "perEdge" else "haplotypes (perBatch)"
        caption = (f"Peak memory (MB) of the algorithms for {mode_text}, small / medium / large datasets: the "
                   "peak minus the memory of the program before the input is read; Ratio = algorithm / static; "
                   "DNF = did not finish within 2 days.")
        lines += [r"\end{tabular}", r"\end{adjustbox}", rf"\caption{{{caption}}}", rf"\label{{tab:memory_{tag}}}",
                  r"\end{table}"]
        out = os.path.join(outdir, f"table_memory_{tag}.tex")
        with open(out, "w") as f:
            f.write("\n".join(lines) + "\n")
        print(f"wrote {out}")


TOOLS = ["bubblegun", "bubblefinder"]
# BENCHMARK files from before 19 Oct 2026 (the first chr22 experiments, older versions of the
# programs: MinI chr22 took 41 h then, 9 h now) are not used
STALE_BEFORE = datetime.datetime(2026, 10, 19).timestamp()


def tool_value(c, tool, mode):
    """the total time of a tool on a dataset (s); DNF if it was not run or exceeds the limit"""
    out = c[2].resolve().parent.parent.parent                       # output/<dataset>
    f = out / tool / MODES[mode][0] / f"{tool}_time.txt"
    if not f.exists():
        return DNF
    t = read_tool_time(f).get("TotalTime")
    return DNF if t is None or t > LIMIT_S else t


def own_value(c, a, mode):
    """(time, starred) of one of our algorithms: the BENCHMARK wall time, else the timed steps"""
    out = c[2].resolve().parent.parent.parent
    f = out / PIPELINE_DIR / MODES[mode][0] / f"{a}_time.txt"
    t = read_pipeline_time(f).get("AlgorithmTime") if f.exists() and f.stat().st_mtime >= STALE_BEFORE else None
    if t is not None:
        return (DNF if t > LIMIT_S else t), False
    v = value(a, c, "TotalTime")
    return v, (v != DNF and not pd.isna(v))


def build_tools(cols_by_mode, outdir):
    """one table per mode: the tools and our algorithms, speedups over BubbleGun"""
    for mode in ("perEdge", "perBatch"):
        cols = cols_by_mode[mode]
        mode_dir, algos, tag = MODES[mode]
        layout = layout_of(cols)

        def row(getter):
            """values and stars of a row, with the group totals"""
            vals, stars, per_group = [], [], {}
            for kind, group, c in layout:
                if kind == "data":
                    v, st = getter(c)
                    per_group.setdefault(group, []).append((v, st))
                    vals.append(v)
                    stars.append(st)
                else:
                    ds = per_group[group]
                    vals.append(group_total([v for v, _ in ds]))
                    stars.append(any(st for _, st in ds))
            return vals, stars
        entries = [(TOOL_LABELS.get(t, t), row(lambda c, t=t: (tool_value(c, t, mode), False)))
                   for t in TOOLS]
        entries += [(name_of(a), row(lambda c, a=a: own_value(c, a, mode))) for a in algos]
        base = entries[0][1][0]
        lines = head_lines(layout)
        lines.append(rf"\texttt{{{entries[0][0]}}} & " + " & ".join(fmt(v, 2) for v in base) + r" \\")
        lines.append(r"\hline")
        for name, (vals, stars) in entries[1:]:
            cells = [fmt(v, 2) + (r"$^*$" if st and v != DNF else "") for v, st in zip(vals, stars)]
            lines.append(rf"\texttt{{{name}}} & " + " & ".join(cells) + r" \\")
            lines.append("Speedup & " + " & ".join(fmt(ratio(b, v), 2) for b, v in zip(base, vals)) + r" \\")
            lines.append(r"\hline")
        mode_text = "edge updates (perEdge)" if mode == "perEdge" else "haplotypes (perBatch)"
        extra = (" The tools did not finish on the HPRC datasets (BubbleGun on chr22 estimated at 91 days, "
                 "BubbleFinder at 42 days)." if mode == "perEdge" else "")
        caption = (f"The external tools and our algorithms for {mode_text}, small / medium / large datasets: "
                   "total time (s), our algorithms as the wall time of the whole run, * = sum of the timed steps; "
                   "speedup = BubbleGun / algorithm; Total = sum over the group; DNF = did not finish within "
                   f"2 days.{extra}")
        lines += [r"\end{tabular}", r"\end{adjustbox}", rf"\caption{{{caption}}}", rf"\label{{tab:tools_{tag}}}",
                  r"\end{table}"]
        out = os.path.join(outdir, f"table_tools_{tag}.tex")
        with open(out, "w") as f:
            f.write("\n".join(lines) + "\n")
        print(f"wrote {out}")


def build(mode, cols, outdir):
    """the long table of a mode: every step of every algorithm"""
    mode_dir, algos, tag = MODES[mode]
    layout = layout_of(cols)
    n_cols = len(layout)
    lines = head_lines(layout)
    for si, (step, suffix) in enumerate(STEPS):
        lines += [rf"\textbf{{{STEP_LABEL[step]}}}" + " &" * n_cols + r" \\", r"\hline"]
        base = row_values(layout, "static", suffix)
        for ai, a in enumerate(algos):
            vals = row_values(layout, a, suffix)
            lines.append(rf"\texttt{{{name_of(a)}}} & " + " & ".join(fmt(v, 4) for v in vals) + r" \\")
            if a != "static":
                lines.append("Speedup & " + " & ".join(fmt(ratio(b, v), 2) for b, v in zip(base, vals)) + r" \\")
            lines.append(r"\hline" if (ai < len(algos) - 1 or si == len(STEPS) - 1) else r"\hline\hline")
    mode_text = "edge updates (perEdge)" if mode == "perEdge" else "haplotypes (perBatch)"
    caption = (f"Time (s) of every step for {mode_text}, small / medium / large datasets; speedup = static / "
               "algorithm; Insert = update of the data structures; Total = sum over the group; DNF = did not "
               "finish within 2 days.")
    lines += [r"\end{tabular}", r"\end{adjustbox}", rf"\caption{{{caption}}}", rf"\label{{tab:full_time_{tag}}}",
              r"\end{table}"]
    out = os.path.join(outdir, f"table_full_time_{tag}.tex")
    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {out}")


def main():
    ap = argparse.ArgumentParser(description="Long time table in three dataset groups")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    cols_by_mode = {}
    for mode in ("perEdge", "perBatch"):
        cols = load(mode)                          # the (large) CSVs are read once per mode
        cols_by_mode[mode] = cols
        build(mode, cols, args.outdir)
        build_steps(mode, cols, args.outdir)
    build_activity(cols_by_mode, args.outdir)
    build_tools(cols_by_mode, args.outdir)
    build_memory(cols_by_mode, args.outdir)


if __name__ == "__main__":
    main()
