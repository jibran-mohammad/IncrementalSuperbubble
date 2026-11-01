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
    table_vertex_processed.tex    explored / affected vertices and what static
    table_edge_processed.tex      processes (graph size summed over the edge
                                  updates / genomes), with ratios; datasets as
                                  columns (csvToLatex.build_activity_totals)
    table_tool_<tool>.tex         one per external tool that has results
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
                                      build_activity_totals, compute_speedup, detect_mode_from_path,
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
    ap.add_argument("--missing-dnf", action="store_true",
                    help="an algorithm without a time file for a dataset is DNF (not run: estimated to "
                         "exceed the limit); with --algorithms, also algorithms without any data")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    mode = detect_mode_from_path(args.csv[0])                      # perEdge / perBatch
    mode_text = "all edge updates (perEdge)" if mode == "perEdge" else "all genomes (perBatch)"
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
        f"algorithm its speedup = static / algorithm. {TIME_NOTE}{DNF_NOTE}{notes}",
        f"tab:paper_components_{tag}")

    # 2. datasets: rows algorithms, columns datasets + Total, total time
    write_table(
        os.path.join(args.outdir, "table_dataset_summary.tex"),
        ["Dataset"] + names + ["Total"],
        algorithm_rows(algos,
                       lambda a, c: total(a, "TotalTime") if c == "Total" else one(a, c, "TotalTime"),
                       names + ["Total"], 4, ratio),
        f"Total time (s), summed over {mode_text}, per dataset and over all datasets; below each "
        f"algorithm its speedup = static / algorithm. {TIME_NOTE}{DNF_NOTE}{notes}",
        f"tab:paper_datasets_{tag}")

    # 3. vertices / edges processed (the transposed activity tables)
    for what, word, fname in (("V", "vertices", "table_vertex_processed.tex"), ("E", "edges", "table_edge_processed.tex")):
        df, cols, algo = build_activity_totals(args.csv, what)
        if algo is None:
            continue
        who = ("the batch algorithm" if algo == "incBSTBatch"
               else "the edge-by-edge algorithms (all give the same counts)")
        generate_latex_table(
            df, cols,
            caption=(f"Total {word} touched, summed over {mode_text}, by {who}: Explored = reached by the "
                     f"searches, Affected = positions reordered; Static per edge / per genome = the graph's "
                     f"{word} summed over the edge updates / genomes, i.e. what static processes in perEdge / "
                     "perBatch mode. Ratios relative to Explored. M = million, B = billion."),
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
            f"incremental DFS and superbubble computation. {TIME_NOTE}{DNF_NOTE}"
            + abbreviation_note([SHORT_NAMES.get(a, a) for a in ds_algos]),
            f"tab:paper_datastructures_{tag}")

    # 5. external tools: tool time vs every algorithm's BENCHMARK wall time
    mode_dir = MODE_DIRS[mode]
    for tool in TOOLS:
        tool_time, ours, parts = {}, {}, {}
        for path, name in zip(args.csv, names):
            out = pathlib.Path(path).resolve().parent.parent.parent      # output/<dataset>
            tf = out / tool / mode_dir / f"{tool}_time.txt"
            if tf.exists():
                tool_time[name] = read_tool_time(tf).get("TotalTime")
                parts[name] = tool_time_parts(tf)
            for a in algos:
                pf = out / PIPELINE_DIR / mode_dir / f"{a}_time.txt"
                if pf.exists():
                    ours[(a, name)] = read_pipeline_time(pf).get("AlgorithmTime")
        if not tool_time:
            continue
        # with --missing-dnf a tool that was not run on a dataset is DNF there
        cols = names if args.missing_dnf else [n for n in names if n in tool_time]
        for n in cols:
            tool_time.setdefault(n, LIMIT_S + 1)
        # the tool table shows the tool, static and the best algorithm of the
        # mode: H1 (incMinIBSTH1) per edge, the batch heuristic per genome
        # (incBSTBatchH1; incBSTBatch while it has no BENCHMARK times yet)
        best = "incMinIBSTH1" if mode == "perEdge" else "incBSTBatchH1"
        if mode != "perEdge" and not any((best, n) in ours for n in cols):
            best = "incBSTBatch"
        # perEdge on HPRC: H1 has no BENCHMARK run on every chromosome; there
        # its time comes from the timed run (sum of its steps), marked *
        fallback = set()
        if mode == "perEdge":
            for n in cols:
                if (best, n) not in ours and finished.get((best, n)) and best in algos:
                    v = one(best, n, "TotalTime")
                    if isinstance(v, (int, float)) and not pd.isna(v):
                        ours[(best, n)] = v
                        fallback.add(n)
        bench_algos = [a for a in ("static", best) if any((a, n) in ours for n in cols)]
        if not bench_algos:
            print(f"[!] {tool}: no BENCHMARK times of our algorithms yet (run pipelineRunner.py / the BENCHMARK runs)")
            continue

        def val(a, c):
            if c == "Total":
                parts = [val(a, n) for n in cols]
                if DNF in parts:
                    return DNF
                return sum(parts) if all(p is not None for p in parts) else None
            v = ours.get((a, c))
            return DNF if v is not None and v > LIMIT_S else v

        tool_total = DNF if any(tool_time[n] is None or tool_time[n] > LIMIT_S for n in cols) \
            else sum(tool_time[n] for n in cols)
        label = TOOL_LABELS.get(tool, tool)
        tool_val = {n: DNF if tool_time[n] is not None and tool_time[n] > LIMIT_S else tool_time[n] for n in cols}
        rows = [[rf"\texttt{{{label}}}"] + [fmt(tool_val[n], 2) for n in cols] + [fmt(tool_total, 2)], None]
        for a in bench_algos:
            vals = [val(a, c) for c in cols + ["Total"]]
            tools = [tool_val[n] for n in cols] + [tool_total]
            star = [("$^*$" if a == best and n in fallback else "") for n in cols] + [""]
            rows.append([rf"\texttt{{{SHORT_NAMES.get(a, a)}}}"] + [fmt(v, 2) + st for v, st in zip(vals, star)])
            rows.append(["Speedup"] + [fmt(ratio(t, v) if v else None, 2) for t, v in zip(tools, vals)])
            rows.append(None)
        fallback_note = (f" $^*$ {SHORT_NAMES.get(best, best)}: no BENCHMARK run; time of the timed run (sum of "
                         f"its steps, without reading the input)." if fallback else "")
        # how much of the tool's time is reading the graph, and the speedup of
        # the best algorithm against the tool's search alone (lower bound),
        # over the datasets with complete (not sampled, not DNF) tool runs
        # and a time of the best algorithm
        ok = [n for n in cols if n in parts and not parts[n][1] and tool_time[n] is not None
              and tool_time[n] <= LIMIT_S and isinstance(val(best, n), (int, float))]
        read_note = ""
        if ok:
            t_all = sum(parts[n][0]["Time"] for n in ok)
            t_read = sum(parts[n][0]["ReadTime"] for n in ok)
            t_find = sum(parts[n][0]["FindTime"] for n in ok)
            t_best = sum(val(best, n) for n in ok)
            if t_all > 0 and t_best > 0:
                read_note = (f" Reading the graph is {100 * t_read / t_all:.0f}\\% of {label}'s time here; "
                             f"counting only its search (no reading), {SHORT_NAMES.get(best, best)} is still "
                             f"{t_find / t_best:.1f}$\\times$ faster (lower bound; total over "
                             f"{'all datasets' if len(ok) == len(cols) else ', '.join(ok)}).")
        write_table(
            os.path.join(args.outdir, f"table_tool_{tool}.tex"),
            ["Dataset"] + cols + ["Total"], rows,
            f"{label} compared with our algorithms ({mode_text}): total time (s); below each algorithm its "
            f"speedup = {label} / algorithm. Our algorithms: wall time of the whole run (reading the edge "
            f"list once, building the base graph, processing every update; BENCHMARK build without timers, "
            f"mean of their runs). {label}: after every update (an edge or a genome) the tool reads the whole "
            f"graph again (it cannot change a graph it has read) and finds all superbubbles from scratch; "
            f"counted are reading and parsing the GFA into its graph and finding the bubbles; not counted "
            f"are program start-up, writing its output, "
            + ("building chains and nesting, " if tool == "bubblegun" else "")
            + "and writing the temporary GFA. Preprocessing (GFA to edge list) is counted on neither side."
            + read_note + fallback_note + DNF_NOTE + abbreviation_note([SHORT_NAMES.get(a, a) for a in bench_algos]),
            f"tab:paper_tool_{tool}_{tag}")


if __name__ == "__main__":
    main()
