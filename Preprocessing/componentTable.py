#!/usr/bin/env python3
"""
componentTable.py - the component summary of the paper: one row per algorithm
(Simple, ImpDFS, ImpDS, ImpSB, with the speedup over Simple), one column per
step, summed over the small datasets. One table per mode.

Times: the sum of the timed steps of every update, from the time files
(<dataset>/timeFile/<mode>/<ID>_time.txt; finished = ends with PeakMemory).

Run from the repository root:
    python3 Preprocessing/componentTable.py --outdir output/paperTables/grouped
writes <outdir>/table_component_summary_<edge|batch>_small.tex and
<outdir>/table_dataset_summary_<edge|batch>_small.tex (total time per dataset),
<outdir>/table_tools_batch_imp.tex (perBatch, the tools and our algorithms; the tools'
reading and search apart, and the speedup over BubbleGun's search alone),
<outdir>/table_memory_batch_imp.tex, table_memory_tools_batch.tex, table_memory_tools_edge.tex
"""
import argparse
import csv
import os
import pathlib
import re
import sys

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "Preprocessing"))
from toolComparisonToLatex import algorithm_time_split, tool_time_split        # noqa: E402
CHROMS = ["11", "17", "18", "22", "X"]
CLASSES = [("small", "small", [f"Human/chromosome_{c}" for c in CHROMS])]
MODES = {"edge": ("perEdge", "edge updates", [("Simple", "static"), ("ImpDFS", "incRecompute"), ("ImpDS", "incMinIBST"),
                                              ("ImpSB", "incMinIBSTH1")]),
         "batch": ("perbatch", "batch updates", [("Simple", "static"), ("ImpDFS", "incRecomputeBatch"),
                                                 ("ImpDS", "incBSTBatch"), ("ImpSB", "incBSTBatchH2")])}
HEAD = ["Insert", "DFS", "Top. Order", "MinI", "MaxO", "SBCompute", "Total"]


def steps(path):
    """the seven step times of a finished run, summed over its updates, else None"""
    if not path.exists():
        return None
    with open(path, "rb") as fh:
        fh.seek(max(0, path.stat().st_size - 4096))
        if b"PeakMemory" not in fh.read():
            return None
    total = [0.0] * 7
    with open(path) as f:
        for line in f:
            p = line.split()
            # update lines: "<number> <genome> <UpdateDS> ... <Total_Time> ..."
            if len(p) >= 10 and p[0].isdigit() and not p[1].replace(".", "").isdigit():
                for i in range(7):
                    total[i] += float(p[2 + i])
    return total


def label(dataset):
    return dataset.split("_")[-1].join(["chr", ""])


def dataset_table(tag, mode_dir, mode_text, algos, cls, datasets, outdir):
    """total time per dataset: one row per algorithm with the speedup over Simple, a Total column"""
    totals = {}
    for name, a in algos:
        vals = [steps(REPO / "output" / d / "timeFile" / mode_dir / f"{a}_time.txt") for d in datasets]
        vals = [v[6] for v in vals]
        totals[name] = vals + [sum(vals)]
    n = len(datasets) + 1
    lines = [r"\begin{table}[ht]", r"\centering", rf"\begin{{tabular}}{{|l|{'|'.join(['c'] * n)}|}}", r"\hline",
             r"\textbf{Dataset} & " + " & ".join(rf"\textbf{{{label(d)}}}" for d in datasets) + r" & \textbf{Total} \\",
             r"\hline"]
    base = totals["Simple"]
    for name, _ in algos:
        lines.append(rf"\texttt{{{name}}} & " + " & ".join(f"{v:.4f}" for v in totals[name]) + r" \\")
        if name != "Simple":
            lines.append("Speedup & " + " & ".join(f"{b / v:.2f}" for b, v in zip(base, totals[name])) + r" \\")
        lines.append(r"\hline")
    lines += [r"\end{tabular}", rf"\caption{{Total time (s) comparison across datasets for {mode_text}}}",
              rf"\label{{tab:total_dataset_summary_{tag}_{cls}}}", r"\end{table}"]
    out = os.path.join(outdir, f"table_dataset_summary_{tag}_{cls}.tex")
    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {out}")


def tool_time(dataset, tool):
    """total time (s) of an external tool, perBatch"""
    for line in open(REPO / "output" / dataset / tool / "perbatch" / f"{tool}_time.txt"):
        k, _, v = line.strip().partition(",")
        if k == "TotalTime":
            return float(v)
    raise SystemExit(f"no TotalTime: {dataset} {tool}")


def tool_split(dataset, tool):
    """{"total", "read", "find"} (s) of an external tool, perBatch: its time, reading the GFA, the search alone"""
    v = tool_time_split(REPO / "output" / dataset / tool / "perbatch" / f"{tool}_time.txt")
    if v is None:
        raise SystemExit(f"no TotalTime: {dataset} {tool}")
    return v


def own_time(dataset, a):
    """(time in s, timed): the BENCHMARK wall time of the whole run, else the sum of the timed steps"""
    f = REPO / "output" / dataset / "incsuperbubble" / "perbatch" / f"{a}_time.txt"
    if f.exists():
        m = re.search(r"AlgorithmTime: ([0-9.]+)", f.read_text())
        if m:
            return float(m.group(1)), False
    return steps(REPO / "output" / dataset / "timeFile" / "perbatch" / f"{a}_time.txt")[6], True


def tools_table(outdir, simple_read=False):
    """perBatch, small / medium / large: the tools and our algorithms, speedups over BubbleGun;
    simple_read: also the row of staticRead (not in the paper yet, local only)"""
    groups = [("Small", [f"Human/chromosome_{c}" for c in CHROMS], True),
              ("Medium", [f"HPRC/chromosome_{c}" for c in CHROMS], True),
              ("Large", ["HPRC2/chromosome_22"], False)]

    def row(get):
        vals = []                                   # (value, timed), with the group totals; None = not run
        for _, ds, total in groups:
            v = [get(d) for d in ds]
            vals += v
            if total:
                vals.append((None if any(x[0] is None for x in v) else sum(x[0] for x in v),
                             any(x[1] for x in v)))
        return vals

    def cells(vals):
        return " & ".join("-" if v is None else f"{v:.2f}" + ("$^*$" if t else "") for v, t in vals)

    def div(b, v):
        return "-" if b is None or not v else f"{b / v:.2f}"

    def speedup(base, vals):
        return "Speedup & " + " & ".join(div(b[0], v[0]) for b, v in zip(base, vals)) + r" \\"

    spec = "|l|" + "||".join("|".join(["c"] * (len(ds) + (1 if t else 0))) for _, ds, t in groups) + "|"
    heads, names = [], []
    for i, (g, ds, t) in enumerate(groups):
        n = len(ds) + (1 if t else 0)
        heads.append(rf"\multicolumn{{{n}}}{{c{'|' if i == len(groups) - 1 else '||'}}}{{\textbf{{{g}}}}}")
        names += ["chr22 v2" if d.startswith("HPRC2") else label(d) for d in ds] + (["Total"] if t else [])
    lines = [r"\begin{table}[ht]", r"\centering", r"\begin{adjustbox}{max width=\linewidth}",
             rf"\begin{{tabular}}{{{spec}}}", r"\hline", " & " + " & ".join(heads) + r" \\", r"\hline",
             r"\textbf{Dataset} & " + " & ".join(rf"\textbf{{{n}}}" for n in names) + r" \\", r"\hline\hline"]
    def split_rows(t):
        """reading the GFA and the search alone below a tool"""
        return [f"{text} & " + cells(row(lambda d: (tool_split(d, t)[k], False))) + r" \\"
                for text, k in (("Reading", "read"), ("Search only", "find"))]

    def search_speedup(vals):
        return "Speedup, search only & " + " & ".join(div(b[0], v[0]) for b, v in zip(gun_find, vals)) + r" \\"

    def static_read(d, k):
        """staticRead (static reading the GFA after every update, as the tools): total, read or find; None = not run"""
        v = algorithm_time_split(REPO / "output" / d / "timeFile" / "perbatch" / "staticRead_time.txt")
        return (None if v is None else v[k]), False

    gun = row(lambda d: (tool_time(d, "bubblegun"), False))
    gun_find = row(lambda d: (tool_split(d, "bubblegun")["find"], False))
    lines += [r"\texttt{BubbleGun} & " + cells(gun) + r" \\"] + split_rows("bubblegun") + [r"\hline"]
    finder = row(lambda d: (tool_time(d, "bubblefinder"), False))
    lines += ([r"\texttt{BubbleFinder} & " + cells(finder) + r" \\"] + split_rows("bubblefinder") +
              [speedup(gun, finder), search_speedup(row(lambda d: (tool_split(d, "bubblefinder")["find"], False))),
               r"\hline"])
    # Simple-Read: static, but after every haplotype it reads the graph from the GFA the tools get
    sread = row(lambda d: static_read(d, "total"))
    lines += [] if not simple_read else [r"\texttt{Simple-Read} & " + cells(sread) + r" \\",
              "Reading & " + cells(row(lambda d: static_read(d, "read"))) + r" \\",
              "Search only & " + cells(row(lambda d: static_read(d, "find"))) + r" \\",
              speedup(gun, sread), search_speedup(row(lambda d: static_read(d, "find"))), r"\hline"]
    for name, a in MODES["batch"][2]:
        vals = row(lambda d: own_time(d, a))
        timed = row(lambda d: (steps(REPO / "output" / d / "timeFile" / "perbatch" / f"{a}_time.txt")[6], True))
        lines += [rf"\texttt{{{name}}} & " + cells(vals) + r" \\", speedup(gun, vals), search_speedup(timed),
                  r"\hline"]
    lines += [r"\end{tabular}", r"\end{adjustbox}",
              r"\caption{The external tools and our algorithms for haplotypes (perBatch), small / medium / large "
              r"datasets: total time (s), our algorithms as the wall time of the whole run, * = sum of the timed "
              r"steps; speedup = BubbleGun / algorithm; Total = sum over the group. Reading / Search only: the "
              r"tool's time split into reading and parsing the GFA (after every haplotype) and the search alone. "
              r"Speedup, search only = BubbleGun's search alone / the sum of the timed steps of our algorithm "
              r"(neither side reads the graph)." +
              (r" Simple-Read = Simple, but after every haplotype it reads the graph from the same GFA as the "
               r"tools into a new adjacency list (Reading) and builds every data structure from scratch (Search "
               r"only; its times are the timed steps, writing the GFA is not counted, as for the tools); - = "
               r"not run." if simple_read else "") + "}",
              r"\label{tab:tools_batch}", r"\end{table}"]
    out = os.path.join(outdir, "table_tools_batch_imp.tex")
    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {out}")


def peak_kb(dataset, mode_dir, a):
    """(peak, baseline) in KB of one of our algorithms, from peak_memory.csv"""
    for row in csv.DictReader(open(REPO / "output" / dataset / "timeFile" / mode_dir / "peak_memory.csv")):
        if row["Algorithm"] == a:
            return float(row["PeakMemoryKB"]), float(row["BaselineMemoryKB"])
    raise SystemExit(f"no memory: {dataset} {a}")


def tool_peak_kb(dataset, mode_dir, tool):
    """whole-process peak memory (KB) of an external tool"""
    for line in open(REPO / "output" / dataset / tool / mode_dir / f"{tool}_time.txt"):
        k, _, v = line.strip().partition(",")
        if k == "PeakMemory(KB)":
            return float(v)
    raise SystemExit(f"no memory: {dataset} {tool}")


def memory_table(out, groups, rows, caption, tag):
    """rows: (name, getter of MB per dataset); the first row is the base of the ratios (no ratio row)"""
    spec = "|l|" + "||".join("|".join(["c"] * len(ds)) for _, ds in groups) + "|"
    all_ds = [d for _, ds in groups for d in ds]
    lines = [r"\begin{table}[ht]", r"\centering", r"\begin{adjustbox}{max width=\linewidth}",
             rf"\begin{{tabular}}{{{spec}}}", r"\hline"]
    heads = [rf"\multicolumn{{{len(ds)}}}{{c{'|' if i == len(groups) - 1 else '||'}}}{{\textbf{{{g}}}}}"
             for i, (g, ds) in enumerate(groups)]
    names = ["chr22 v2" if d.startswith("HPRC2") else label(d) for d in all_ds]
    lines += [" & " + " & ".join(heads) + r" \\", r"\hline",
              r"\textbf{Dataset} & " + " & ".join(rf"\textbf{{{n}}}" for n in names) + r" \\", r"\hline\hline"]
    base = None
    for name, get in rows:
        vals = [get(d) for d in all_ds]
        lines.append(rf"\texttt{{{name}}} & " + " & ".join(f"{v:.1f}" for v in vals) + r" \\")
        if base is None:
            base = vals
        else:
            lines.append("Ratio & " + " & ".join(f"{v / b:.2f}" for v, b in zip(vals, base)) + r" \\")
        lines.append(r"\hline")
    lines += [r"\end{tabular}", r"\end{adjustbox}", rf"\caption{{{caption}}}", rf"\label{{{tag}}}", r"\end{table}"]
    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {out}")


def memory_tables(outdir):
    small = ("Small", [f"Human/chromosome_{c}" for c in CHROMS])
    all_groups = [small, ("Medium", [f"HPRC/chromosome_{c}" for c in CHROMS]), ("Large", ["HPRC2/chromosome_22"])]

    def net(mode_dir, a):
        def get(d):
            peak, base = peak_kb(d, mode_dir, a)
            return (peak - base) / 1024.0
        return get

    def whole(mode_dir, a):
        return lambda d: peak_kb(d, mode_dir, a)[0] / 1024.0

    def tool(mode_dir, t):
        return lambda d: tool_peak_kb(d, mode_dir, t) / 1024.0

    batch, edge = MODES["batch"], MODES["edge"]
    memory_table(os.path.join(outdir, "table_memory_batch_imp.tex"), all_groups,
                 [(n, net(batch[0], a)) for n, a in batch[2]],
                 "Peak memory (MB) of the algorithms for haplotypes (perBatch), small / medium / large datasets: the "
                 "peak minus the memory of the program before the input is read; Ratio = algorithm / Simple.",
                 "tab:memory_batch")
    for tag, (mode_dir, _, algos), groups, text in (("batch", batch, all_groups, "haplotypes (perBatch), small / medium / large datasets"),
                                                    ("edge", edge, [small], "edge updates (perEdge), small datasets")):
        rows = [("BubbleGun", tool(mode_dir, "bubblegun")), ("BubbleFinder", tool(mode_dir, "bubblefinder"))]
        rows += [(n, whole(mode_dir, a)) for n, a in algos]
        memory_table(os.path.join(outdir, f"table_memory_tools_{tag}.tex"), groups, rows,
                     f"Peak memory (MB) of the external tools and our algorithms for {text}: the whole-process peak; "
                     "Ratio = tool or algorithm / BubbleGun.", f"tab:memory_tools_{tag}")


def main():
    ap = argparse.ArgumentParser(description="Component summary per mode and class")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--simple-read", action="store_true",
                    help="also the Simple-Read row (staticRead) in the tool table; not in the paper yet")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    tools_table(args.outdir, args.simple_read)
    memory_tables(args.outdir)
    for tag, (mode_dir, mode_text, algos) in MODES.items():
        for cls, cls_text, datasets in CLASSES:
            per_algo = {name: [0.0] * 7 for name, _ in algos}
            used = []
            for d in datasets:
                runs = {name: steps(REPO / "output" / d / "timeFile" / mode_dir / f"{a}_time.txt") for name, a in algos}
                if any(r is None for r in runs.values()):
                    print(f"  {tag} {cls}: {d} left out (not all finished: "
                          f"{[n for n, r in runs.items() if r is None]})")
                    continue
                used.append(d)
                for name, r in runs.items():
                    per_algo[name] = [x + y for x, y in zip(per_algo[name], r)]
            if len(used) == len(datasets):
                dataset_table(tag, mode_dir, mode_text, algos, cls, datasets, args.outdir)
            if not used:
                print(f"  {tag} {cls}: no table, no dataset where all four finished")
                continue
            lines = [r"\begin{table}[ht]", r"\centering", r"\begin{adjustbox}{max width=\linewidth}",
                     r"\begin{tabular}{|l|c|c|c|c|c|c|c|}", r"\hline",
                     r"\textbf{Component} & " + " & ".join(rf"\textbf{{{h}}}" for h in HEAD) + r" \\", r"\hline"]
            base = per_algo["Simple"]
            for name, _ in algos:
                lines.append(rf"\texttt{{{name}}} & " + " & ".join(f"{v:.4f}" for v in per_algo[name]) + r" \\")
                if name != "Simple":
                    lines.append("Speedup & " + " & ".join(f"{b / v:.2f}" if v > 0 else "-"
                                                           for b, v in zip(base, per_algo[name])) + r" \\")
                lines.append(r"\hline")
            on = "" if len(used) == len(datasets) else " (" + ", ".join(label(d) for d in used) + " only)"
            caption = f"Total time (s) for each component for {mode_text}, {cls_text} datasets{on}"
            lines += [r"\end{tabular}", r"\end{adjustbox}", rf"\caption{{{caption}}}",
                      rf"\label{{tab:total_component_summary_{tag}_{cls}}}", r"\end{table}"]
            out = os.path.join(args.outdir, f"table_component_summary_{tag}_{cls}.tex")
            with open(out, "w") as f:
                f.write("\n".join(lines) + "\n")
            print(f"wrote {out}")


if __name__ == "__main__":
    main()
