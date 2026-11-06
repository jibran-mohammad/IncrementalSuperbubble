"""
paperTables.py - the tables of the paper, kept apart from the detailed
tables of csvToLatex.py (output/latexTables*).

For one mode (the CSVs of parse_to_csv.py, one per dataset) it writes to
--outdir:

    table_component_summary.tex   rows: algorithms (static first, below every
                                  other algorithm its speedup over static);
                                  columns: the steps (UpdateDS ... SBCompute)
                                  and Total; times summed over all datasets
    table_dataset_summary.tex     rows: the same; columns: the datasets and
                                  Total; total time per dataset
    table_full_time.tex           rows: one block per step (UpdateDS ...
                                  SBCompute, Total), in each the algorithms
                                  with their speedup over static; columns: the
                                  datasets and Total
    table_step_<step>.tex         the same for one step (UpdateDS, DFS,
                                  TopoOrder, MinI, MaxO, SBCompute)
    table_memory.tex              rows: algorithms (below each its memory
                                  relative to static); columns: the datasets;
                                  peak memory minus the program's own memory
                                  before reading the input (peak_memory.csv)
    table_size.tex                (with --size) genomes and graph size (base
                                  and final graph) per dataset
    table_vertex_processed.tex    explored / affected vertices and what static
    table_edge_processed.tex      processes (graph size summed over the edge
                                  updates / genomes), with ratios; datasets as
                                  columns (csvToLatex.build_activity_totals)
    table_tools.tex               the external tools that have results
                                  (output/<dataset>/<tool>/<mode>/<tool>_time.txt):
                                  the tool's time, then static and the best
                                  algorithm of the mode (perEdge: H1, perBatch:
                                  the batch algorithm), each with its speedup =
                                  tool / algorithm; datasets as
                                  columns. Both sides are wall times of the
                                  whole program: the algorithms from their
                                  BENCHMARK runs (output/<dataset>/incsuperbubble/
                                  <mode>/<ID>_time.txt, AlgorithmTime), the
                                  tool from its own run (TotalTime).

Times in the first two tables are the per-update timers of the normal build
(as in csvToLatex.py), the tool table uses the BENCHMARK wall times; the
captions say which.

Usage (from the repository root):
    python3 Preprocessing/paperTables.py --csv <CSV> ... --outdir output/paperTables/perEdge
"""
import argparse
import os
import pathlib
import sys

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))
from Preprocessing.csvToLatex import (ALGORITHMS, SHORT_NAMES, TIME_NOTE, abbreviation_note,
                                      build_activity_totals, build_size_haplotype_table, compute_speedup, detect_mode_from_path,
                                      generate_latex_table, infer_dataset_name_from_path,
                                      normalize_dataframe, sum_metric)
from Preprocessing.toolComparisonToLatex import TOOL_LABELS, read_pipeline_time, read_tool_time
from pipelineRunner import MODE_DIRS, TOOL_NAME as PIPELINE_DIR

STEPS = [("UpdateDS", "UpdateDSTime"), ("DFS", "DFSTime"), ("TopoOrder", "TopologicalOrderTime"),
         ("MinI", "MinITime"), ("MaxO", "MaxOTime"), ("SBCompute", "SBComputeTime"),
         ("Total", "TotalTime")]
TOOLS = ["bubblegun", "bubblefinder"]

# Scalability rule: a run (ours, static or a tool) that needs more than
# LIMIT_DAYS days - measured, stopped at the limit, or estimated - is shown
# as DNF ("did not finish within LIMIT_DAYS days") and gets no speedup.
LIMIT_DAYS = 2
LIMIT_S = LIMIT_DAYS * 86400
DNF = "DNF"
DNF_NOTE = (f" DNF = did not finish within {LIMIT_DAYS} days (not scalable): stopped at the limit, or "
            f"its measured or estimated time exceeds {LIMIT_DAYS} days.")


def run_finished(csv_path, algo, missing_dnf=False):
    """False if the algorithm's time file next to the CSV ends without its
    final PeakMemory line, i.e. the run was stopped (time limit) or crashed;
    with missing_dnf also if there is no time file (not run: estimated to
    exceed the limit)."""
    f = pathlib.Path(csv_path).parent / f"{algo}_time.txt"
    if not f.exists():
        return not missing_dnf
    with open(f, "rb") as fh:
        fh.seek(max(0, f.stat().st_size - 4096))
        return b"PeakMemory" in fh.read()


def fmt(x, digits):
    if x is None or (isinstance(x, float) and pd.isna(x)):
        return "-"
    if x == DNF:
        return DNF
    return f"{x:.{digits}f}"


def ratio(base, value):
    """speedup base / value; none if either side did not finish"""
    if base in (None, DNF) or value in (None, DNF):
        return None
    return compute_speedup(base, value)


def tool_time_parts(path):
    """Sums of Time / ReadTime / FindTime over the rows (base graph and every
    update) of a <tool>_time.txt, and whether TotalTime is an estimate
    (sampled perEdge run: the rows are only the sampled updates)."""
    sums, idx, estimated = {"Time": 0.0, "ReadTime": 0.0, "FindTime": 0.0}, None, False
    for line in open(path):
        f = line.strip().split(",")
        if idx is None and "ReadTime" in f:
            idx = {k: f.index(k) for k in sums if k in f}
            continue
        if f[0] == "TotalTimeEstimated":
            estimated = True
        if idx and (f[0].isdigit() or f[0] == "Base"):
            for k, i in idx.items():
                try:
                    sums[k] += float(f[i])
                except (ValueError, IndexError):
                    pass
    return sums, estimated


def write_table(path, header, rows, caption, label):
    """header: column titles; rows: lists of cells, None = \\hline"""
    lines = [r"\begin{table}[ht]", r"\centering",
             r"\begin{tabular}{|l|" + "|".join(["c"] * (len(header) - 1)) + "|}", r"\hline",
             " & ".join(rf"\textbf{{{h}}}" for h in header) + r" \\", r"\hline"]
    for r in rows:
        lines.append(r"\hline" if r is None else " & ".join(r) + r" \\")
    lines += [r"\end{tabular}", rf"\caption{{{caption}}}", rf"\label{{{label}}}", r"\end{table}"]
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {path}")


def base_graph_size(csv_path):
    """(vertices, edges) of the base graph without the artificial source and
    sink: the line after the graph name in input/<dataset>/Graphs/
    numeric_edgelist.txt ("<edges> <vertices>"); without the edge list the
    first row (step 0, the base graph) of a history summary of the dataset
    (output/<dataset>/history/*_history_summary.tsv: new_edges,
    new_vertices); None if neither is found"""
    out = pathlib.Path(csv_path).resolve().parent.parent.parent        # output/<dataset>
    rel = out.relative_to(REPO / "output")
    f = REPO / "input" / rel / "Graphs" / "numeric_edgelist.txt"
    if not f.exists():
        for h in sorted((out / "history").glob("*_history_summary.tsv")):
            with open(h) as fh:
                head = fh.readline().rstrip("\n").split("\t")
                row = dict(zip(head, fh.readline().rstrip("\n").split("\t")))
            if row.get("step") == "0" and "new_edges" in row:
                return int(row["new_vertices"]), int(row["new_edges"])
        return None
    with open(f) as fh:
        line = fh.readline()
        while line and (not line.strip() or line.startswith("#")):
            line = fh.readline()
        e, v = fh.readline().split()[:2]
    return int(v), int(e)


def history_graph_size(csv_path):
    """(vertices, edges) of the final graph without the artificial source and
    sink from a history summary of the dataset: new_vertices / new_edges
    summed over all steps (the base graph is step 0); None if there is none"""
    out = pathlib.Path(csv_path).resolve().parent.parent.parent        # output/<dataset>
    for h in sorted((out / "history").glob("*_history_summary.tsv")):
        v = e = 0
        with open(h) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            if "new_edges" not in head:
                continue
            ce, cv = head.index("new_edges"), head.index("new_vertices")
            for line in fh:
                p = line.rstrip("\n").split("\t")
                e += int(p[ce]); v += int(p[cv])
        return v, e
    return None


def algorithm_rows(algos, value_of, columns, digits, ratio_fn, base_algo="static"):
    """base_algo (static) first, then every other algorithm with a Speedup row
    below it. value_of(algo, column) -> number; ratio(base, value) -> speedup."""
    rows = []
    for a in algos:
        vals = [value_of(a, c) for c in columns]
        rows.append([rf"\texttt{{{SHORT_NAMES.get(a, a)}}}"] + [fmt(v, digits) for v in vals])
        if a != base_algo:
            base = [value_of(base_algo, c) for c in columns]
            rows.append(["Speedup"] + [fmt(ratio_fn(b, v), 2) for b, v in zip(base, vals)])
        rows.append(None)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", nargs="+", required=True, help="haplotype_build_comparison.csv, one per dataset")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tag", help="label suffix, e.g. perbatch_hprc (default: the mode)")
    ap.add_argument("--algorithms", nargs="+",
                    help="only these algorithm IDs (default: every algorithm with data), "
                         "e.g. static incBSTBatch for all HPRC chromosomes")
    ap.add_argument("--size", action="store_true",
                    help="also write table_size.tex (haplotypes and graph size per dataset; the same in both "
                         "modes, so one mode is enough)")
    ap.add_argument("--tool-note", default="",
                    help="sentence added to the perEdge tool-table captions after the BENCHMARK / timed ratios")
    ap.add_argument("--dnf-note", default="",
                    help="sentence added after the DNF note of every caption, e.g. how the runs that were "
                         "not started were estimated")
    ap.add_argument("--missing-dnf", action="store_true",
                    help="an algorithm without a time file for a dataset is DNF (not run: estimated to "
                         "exceed the limit); with --algorithms, also algorithms without any data")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    dnf_note = DNF_NOTE + (" " + args.dnf_note if args.dnf_note else "")
    mode = detect_mode_from_path(args.csv[0])                      # perEdge / perBatch
    mode_text = "all edge updates (perEdge)" if mode == "perEdge" else "all haplotypes (perBatch)"
    tag = args.tag or mode.lower()
    names = [infer_dataset_name_from_path(p) for p in args.csv]
    frames = [normalize_dataframe(pd.read_csv(p, low_memory=False), p) for p in args.csv]
    algos = [a for a in ALGORITHMS if any(not pd.isna(sum_metric(df, a, "TotalTime")) for df in frames)]
    if args.algorithms:
        algos = [a for a in (ALGORITHMS if args.missing_dnf else algos) if a in args.algorithms]
    notes = abbreviation_note([SHORT_NAMES.get(a, a) for a in algos])

    finished = {(a, n): run_finished(p, a, args.missing_dnf) for a in algos for p, n in zip(args.csv, names)}

    def one(a, n, suffix):
        """a step's time of one dataset; DNF if the run did not finish or took
        longer than the limit (the limit applies to a run, not to a sum over
        datasets)"""
        if not finished[(a, n)]:
            return DNF
        run = sum_metric(by_name[n], a, "TotalTime")
        if not pd.isna(run) and run > LIMIT_S:
            return DNF
        return sum_metric(by_name[n], a, suffix)

    def total(a, suffix, frame_list=None):
        vals = [one(a, n, suffix) for n in names]
        if DNF in vals:
            return DNF
        vals = [v for v in vals if not pd.isna(v)]
        return sum(vals) if vals else float("nan")

    by_name = dict(zip(names, frames))
    # 1. components: rows algorithms, columns steps, summed over all datasets
    write_table(
        os.path.join(args.outdir, "table_component_summary.tex"),
        ["Component"] + [s for s, _ in STEPS],
        algorithm_rows(algos, lambda a, c: total(a, dict(STEPS)[c]), [s for s, _ in STEPS], 4, ratio),
        f"Time (s) of every step, summed over {mode_text} of all {len(names)} datasets; below each "
        f"algorithm its speedup = static / algorithm. {TIME_NOTE}{dnf_note}{notes}",
        f"tab:paper_components_{tag}")

    # 2. datasets: rows algorithms, columns datasets + Total, total time
    write_table(
        os.path.join(args.outdir, "table_dataset_summary.tex"),
        ["Dataset"] + names + ["Total"],
        algorithm_rows(algos,
                       lambda a, c: total(a, "TotalTime") if c == "Total" else one(a, c, "TotalTime"),
                       names + ["Total"], 4, ratio),
        f"Total time (s), summed over {mode_text}, per dataset and over all datasets; below each "
        f"algorithm its speedup = static / algorithm. {TIME_NOTE}{dnf_note}{notes}",
        f"tab:paper_datasets_{tag}")

    # 2b. every step: one block per step in one table (table_full_time), and
    # one table per step; rows algorithms, columns datasets + Total
    def step_value(suffix):
        return lambda a, c: total(a, suffix) if c == "Total" else one(a, c, suffix)
    def full_rows(alist):
        rows = []
        for step, suffix in STEPS:
            rows.append([rf"\textbf{{{step}}}"] + [""] * (len(names) + 1))
            rows.append(None)
            rows += algorithm_rows(alist, step_value(suffix), names + ["Total"], 4, ratio)
        return rows
    # perEdge: the algorithms (without NBST and BST) in table_full_time, the
    # data structures for minI / maxO (static, NBST, BST, MinI) in
    # table_full_time_trees, next to it
    trees = ["static", "incNormalBST", "incBST", "incMinIBST"]
    split = mode == "perEdge" and any(a in algos for a in ("incNormalBST", "incBST"))
    full = full_rows([a for a in algos if not (split and a in ("incNormalBST", "incBST"))])
    for step, suffix in STEPS:
        if step != "Total":
            write_table(
                os.path.join(args.outdir, f"table_step_{step}.tex"),
                ["Dataset"] + names + ["Total"],
                algorithm_rows(algos, step_value(suffix), names + ["Total"], 4, ratio),
                f"Time (s) of the step {step}, summed over {mode_text}, per dataset and over all datasets; "
                f"below each algorithm its speedup = static / algorithm. {TIME_NOTE}{dnf_note}{notes}",
                f"tab:paper_step_{step.lower()}_{tag}")
    shown = [a for a in algos if not (split and a in ("incNormalBST", "incBST"))]
    write_table(
        os.path.join(args.outdir, "table_full_time.tex"),
        ["Dataset"] + names + ["Total"], full,
        f"Time (s) of every step, summed over {mode_text}, per dataset and over all datasets; below each "
        f"algorithm its speedup = static / algorithm."
        + (" The data structures for minI and maxO (NBST, BST) are compared in the next table." if split else "")
        + f" {TIME_NOTE}{dnf_note}" + abbreviation_note([SHORT_NAMES.get(a, a) for a in shown]),
        f"tab:paper_full_time_{tag}")
    if split:
        talgos = [a for a in trees if a in algos]
        write_table(
            os.path.join(args.outdir, "table_full_time_trees.tex"),
            ["Dataset"] + names + ["Total"], full_rows(talgos),
            f"The data structures for minI and maxO: time (s) of every step, summed over {mode_text}, per "
            f"dataset and over all datasets; below each its speedup = static / it. NBST: normal balanced "
            f"binary search trees keyed by position (every vertex that moves is re-keyed); BST: modified trees "
            f"that compare by the current positions, for minI and maxO; MinI: modified tree for minI, maxO by "
            f"comparison. All use the same incremental DFS and superbubble computation. {TIME_NOTE}{dnf_note}"
            + abbreviation_note([SHORT_NAMES.get(a, a) for a in talgos]),
            f"tab:paper_full_time_trees_{tag}")

    # 2c. peak memory per dataset: the peak (VmHWM) minus the program's own
    # memory before the input is read, i.e. the graph and the data structures
    memory = {}
    for path, name in zip(args.csv, names):
        f = pathlib.Path(path).parent / "peak_memory.csv"
        if f.exists():
            m = pd.read_csv(f)
            for _, r in m.iterrows():
                base = r["BaselineMemoryKB"] if "BaselineMemoryKB" in m.columns and not pd.isna(r["BaselineMemoryKB"]) else 0
                memory[(r["Algorithm"], name)] = (r["PeakMemoryKB"] - base) / 1024.0

    def mem_value(a, c):
        if not finished[(a, c)]:
            return DNF
        return memory.get((a, c), float("nan"))

    def mem_ratio(base, value):
        return None if base in (None, DNF) or value in (None, DNF) or pd.isna(base) or pd.isna(value) \
            else value / base
    if memory:
        rows = algorithm_rows(algos, mem_value, names, 1, mem_ratio)
        rows = [["Ratio" if r and r[0] == "Speedup" else c for c in r][:1] + r[1:] if r else r for r in rows]
        write_table(
            os.path.join(args.outdir, "table_memory.tex"),
            ["Dataset"] + names, rows,
            f"Peak memory (MB) of each algorithm ({mode_text}): the peak resident memory of the program "
            f"(VmHWM) minus its own memory measured before the input is read (about 3.5 MB), i.e. the "
            f"graph and the algorithm's data structures; below each algorithm its memory relative to "
            f"static (Ratio = algorithm / static).{dnf_note}{notes}",
            f"tab:paper_memory_{tag}")

    # 2d. dataset size (the same in both modes: written with --size only)
    if args.size:
        size = build_size_haplotype_table(args.csv)
        base = {name: base_graph_size(path) for path, name in zip(args.csv, names)}
        # without RealNodes in dataset_info.csv (no edge list here) the final
        # size there includes the artificial source / sink: take the history
        real = {}
        for path, name in zip(args.csv, names):
            info = pathlib.Path(path).parent / "dataset_info.csv"
            if info.exists() and pd.isna(pd.read_csv(info).iloc[0].get("RealNodes")):
                real[name] = history_graph_size(path)
        rows = []
        for _, r in size.iterrows():
            n = r["Dataset"]
            if real.get(n):
                r = r.copy()
                r["TotalNodes"], r["TotalEdges"] = real[n]
            bv, be = base.get(n) or ("-", "-")
            rows.append([n, f"{int(r['Haplotypes']):,}".replace(",", "{,}"),
                         f"{bv:,}".replace(",", "{,}") if bv != "-" else "-",
                         f"{be:,}".replace(",", "{,}") if be != "-" else "-",
                         f"{int(r['TotalNodes']):,}".replace(",", "{,}"),
                         f"{int(r['TotalEdges']):,}".replace(",", "{,}")])
        write_table(os.path.join(args.outdir, "table_size.tex"),
                    ["Dataset", "Haplotypes", "Base vertices", "Base edges", "Final vertices", "Final edges"], rows,
                    "Datasets: haplotypes added to the base graph, and the size of the base graph and of the final "
                    "graph (vertices and edges, without the artificial source and sink).",
                    f"tab:paper_size_{tag}")

    # 3. vertices / edges processed (the transposed activity tables)
    for what, word, fname in (("V", "vertices", "table_vertex_processed.tex"), ("E", "edges", "table_edge_processed.tex")):
        df, cols, algo = build_activity_totals(args.csv, what)
        if algo is None:
            continue
        who = ("the batch algorithm" if algo == "incBSTBatch"
               else "the edge-by-edge algorithms (all give the same counts)")
        # only the static row of this mode (per edge in perEdge, per genome in perBatch)
        other = "per haplotype" if mode == "perEdge" else "per edge"
        df = df[~df["Dataset"].str.contains(other)].reset_index(drop=True)
        own = "edge update" if mode == "perEdge" else "haplotype"
        generate_latex_table(
            df, cols,
            caption=(f"Total {word} touched, summed over {mode_text}, by {who}: Explored = reached by the "
                     f"searches, Affected = positions reordered; Static per {own.split()[0]} = the graph's "
                     f"{word} summed over the {own}s, i.e. what static processes in {mode}. Ratios relative "
                     "to Explored. M = million, B = billion."),
            label=f"tab:paper_{word}_processed_{tag}",
            out_tex=os.path.join(args.outdir, fname))

    # 4. the data structures for minI / maxO (perEdge): normal BST (NBST),
    # modified BST for both (BST), modified BST for minI and maxO by
    # comparison (MinI); speedup relative to NBST
    ds_algos = ["incNormalBST", "incBST", "incMinIBST"]
    if mode == "perEdge" and all(a in algos for a in ds_algos) and \
            any(total(a, "TotalTime") != DNF for a in ds_algos):           # not only DNF (HPRC)
        ds_cols = [("MinI", "MinITime"), ("MaxO", "MaxOTime"), ("MinI+MaxO", None), ("Total", "TotalTime")]

        def ds_value(a, c):
            if c == "MinI+MaxO":
                mi, mo = total(a, "MinITime"), total(a, "MaxOTime")
                return DNF if DNF in (mi, mo) else mi + mo
            return total(a, dict(ds_cols)[c])
        write_table(
            os.path.join(args.outdir, "table_datastructures.tex"),
            ["Data structure"] + [c for c, _ in ds_cols],
            algorithm_rows(ds_algos, ds_value, [c for c, _ in ds_cols], 4, ratio, base_algo="incNormalBST"),
            f"The data structures for minI and maxO: time (s) of maintaining minI and maxO and the total "
            f"time, summed over {mode_text} of all {len(names)} datasets; below each its speedup = "
            f"NBST / algorithm. NBST: normal balanced binary search trees keyed by position (every "
            f"vertex that moves is re-keyed); BST: modified trees that compare by the current positions, "
            f"for minI and maxO; MinI: modified tree for minI, maxO by comparison. All three use the same "
            f"incremental DFS and superbubble computation. {TIME_NOTE}{dnf_note}"
            + abbreviation_note([SHORT_NAMES.get(a, a) for a in ds_algos]),
            f"tab:paper_datastructures_{tag}")

    # 5. external tools: one table, BubbleGun as the reference: BubbleGun,
    # then BubbleFinder, static and the best algorithm of the mode, each with
    # its speedup = BubbleGun / it. Our algorithms: BENCHMARK wall times.
    mode_dir = MODE_DIRS[mode]
    tool_time, parts, ours = {t: {} for t in TOOLS}, {t: {} for t in TOOLS}, {}
    for path, name in zip(args.csv, names):
        out = pathlib.Path(path).resolve().parent.parent.parent      # output/<dataset>
        for tool in TOOLS:
            tf = out / tool / mode_dir / f"{tool}_time.txt"
            if tf.exists():
                tool_time[tool][name] = read_tool_time(tf).get("TotalTime")
                parts[tool][name] = tool_time_parts(tf)
        for a in algos:
            pf = out / PIPELINE_DIR / mode_dir / f"{a}_time.txt"
            if pf.exists():
                ours[(a, name)] = read_pipeline_time(pf).get("AlgorithmTime")
    tools = [t for t in TOOLS if tool_time[t]]
    if tools:
        # with --missing-dnf a tool that was not run on a dataset is DNF there
        cols = names if args.missing_dnf else [n for n in names if any(n in tool_time[t] for t in tools)]
        for t in tools:
            for n in cols:
                tool_time[t].setdefault(n, LIMIT_S + 1)
        # the best algorithm of the mode: H1 (incMinIBSTH1) per edge, the batch
        # heuristic per haplotype (incBSTBatchH1; incBSTBatch while it has no
        # BENCHMARK times yet)
        best = "incMinIBSTH1" if mode == "perEdge" else "incBSTBatchH1"
        if mode != "perEdge" and not any((best, n) in ours for n in cols):
            best = "incBSTBatch"
        # perEdge on HPRC: H1 without a BENCHMARK run on a dataset: the time of
        # the timed run (sum of its steps), marked *
        fallback, both = set(), []
        if mode == "perEdge":
            for n in cols:
                timed = one(best, n, "TotalTime") if best in algos and finished.get((best, n)) else None
                if (best, n) in ours and isinstance(timed, (int, float)) and not pd.isna(timed):
                    both.append((n, timed, ours[(best, n)]))
                if (best, n) not in ours and finished.get((best, n)) and best in algos:
                    v = one(best, n, "TotalTime")
                    if isinstance(v, (int, float)) and not pd.isna(v):
                        ours[(best, n)] = v
                        fallback.add(n)
        bench_algos = [a for a in ("static", best) if any((a, n) in ours for n in cols)]

        def ours_val(a, c):
            if c == "Total":
                ps = [ours_val(a, n) for n in cols]
                if DNF in ps:
                    return DNF
                return sum(ps) if all(p is not None for p in ps) else None
            v = ours.get((a, c))
            return DNF if v is not None and v > LIMIT_S else v

        def tool_val(t, c):
            if c == "Total":
                ps = [tool_val(t, n) for n in cols]
                if DNF in ps:
                    return DNF
                return sum(ps) if all(p is not None for p in ps) else None
            v = tool_time[t].get(c)
            return DNF if v is not None and v > LIMIT_S else v

        base = tools[0]                                       # BubbleGun (or the only tool there is)
        base_label = TOOL_LABELS.get(base, base)
        ref = [tool_val(base, c) for c in cols + ["Total"]]
        rows = [[rf"\texttt{{{base_label}}}"] + [fmt(v, 2) for v in ref], None]
        entries = [(TOOL_LABELS.get(t, t), [tool_val(t, c) for c in cols + ["Total"]], None) for t in tools[1:]]
        entries += [(SHORT_NAMES.get(a, a), [ours_val(a, c) for c in cols + ["Total"]], a) for a in bench_algos]
        for name, vals, a in entries:
            star = [("$^*$" if a == best and n in fallback else "") for n in cols] + [""]
            rows.append([rf"\texttt{{{name}}}"] + [fmt(v, 2) + st for v, st in zip(vals, star)])
            rows.append(["Speedup"] + [fmt(ratio(r, v) if v else None, 2) for r, v in zip(ref, vals)])
            rows.append(None)

        fallback_note = (f" $^*$ {SHORT_NAMES.get(best, best)}: no BENCHMARK run; time of the timed build (the "
                         f"sum of its timed steps, without reading the input; the timed build also writes a "
                         f"time line per update)." if fallback else "")
        # BENCHMARK (wall time of the whole program) vs the timed steps of
        # the step tables: always higher; the measured ratio per dataset
        bench = [(n, t, b) for n, t, b in both if n not in fallback]
        if mode == "perEdge" and bench:
            fallback_note += (f" The time and step tables sum only the timed steps of {SHORT_NAMES.get(best, best)} "
                              f"(without reading the input and the work between the steps), so the wall time "
                              f"here is higher: " + ", ".join(f"{n} {b / t:.2f}$\\times$" for n, t, b in bench)
                              + "." + (" " + args.tool_note if args.tool_note else ""))
        # per tool: how much of its time is reading the graph, and the speedup
        # of the best algorithm against its search alone (lower bound), over
        # the datasets with complete (not sampled, not DNF) tool runs
        read_note = ""
        for t in tools:
            label = TOOL_LABELS.get(t, t)
            ok = [n for n in cols if n in parts[t] and not parts[t][n][1] and tool_time[t][n] is not None
                  and tool_time[t][n] <= LIMIT_S and isinstance(ours_val(best, n), (int, float))]
            if not ok:
                continue
            t_all = sum(parts[t][n][0]["Time"] for n in ok)
            t_read = sum(parts[t][n][0]["ReadTime"] for n in ok)
            t_find = sum(parts[t][n][0]["FindTime"] for n in ok)
            t_best = sum(ours_val(best, n) for n in ok)
            if t_all > 0 and t_best > 0:
                read_note += (f" Reading the graph is {100 * t_read / t_all:.0f}\\% of {label}'s time here; "
                              f"counting only its search (no reading), {SHORT_NAMES.get(best, best)} is still "
                              f"{t_find / t_best:.1f}$\\times$ faster than {label} (lower bound; total over "
                              f"{'all datasets' if len(ok) == len(cols) else ', '.join(ok)}).")
        labels = [TOOL_LABELS.get(t, t) for t in tools]
        write_table(
            os.path.join(args.outdir, "table_tools.tex"),
            ["Dataset"] + cols + ["Total"], rows,
            f"The external tools compared with our algorithms ({mode_text}): total time (s); below each its "
            f"speedup = {base_label} / it. Our algorithms: wall time of the whole run (reading the edge list "
            f"once, building the base graph, processing every update; BENCHMARK build without timers, mean of "
            f"their runs). {' and '.join(labels)}: after every update (an edge or a haplotype) the tool reads "
            f"the whole graph again (it cannot change a graph it has read) and finds all superbubbles from "
            f"scratch; counted are reading and parsing the GFA into its graph and finding the bubbles; not "
            f"counted are program start-up, writing its output, building chains and nesting (BubbleGun), and "
            f"writing the temporary GFA. Preprocessing (GFA to edge list) is counted on neither side."
            + read_note + fallback_note + dnf_note
            + abbreviation_note([SHORT_NAMES.get(a, a) for a in bench_algos]),
            f"tab:paper_tools_{tag}")

if __name__ == "__main__":
    main()
