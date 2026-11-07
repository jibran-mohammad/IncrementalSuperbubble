"""
smallGraphsTest.py

Reproduces the experiments of the paper on the five small human graphs
(chromosomes 11, 17, 18, 22 and X) in one command: our algorithms per edge
and per batch, BubbleGun and BubbleFinder, and tables of the results.

=====================================================================
HOW TO RUN
=====================================================================

Run from the repository root, with the virtual environment active:

    python3 -m venv venv && source venv/bin/activate
    pip install -r requirements.txt

    # all five graphs, per batch tools (about 40 minutes)
    python3 scripts/smallGraphsTest.py

    # one graph only (chr22, the smallest: a few minutes)
    python3 scripts/smallGraphsTest.py -c 22

    # also BubbleGun, BubbleFinder and Simple-Read after every edge: this takes
    # HOURS (about 8 hours for the five graphs: 4.5 for the tools, 3.5 for
    # Simple-Read), see the note below
    python3 scripts/smallGraphsTest.py --tools-per-edge

Options:
    -c / --chromosomes  11 17 18 22 X               (default: all five)
    --runs N            runs of each of our algorithms, averaged (default: 1;
                        the paper used 30)
    --tools-per-edge    also run BubbleGun, BubbleFinder and Simple-Read per
                        edge (hours)
    --no-tools          only our algorithms
    --overwrite         run even if results already exist in output/ (they
                        are replaced)

Per edge, BubbleGun and BubbleFinder cannot update their result: after every
inserted edge they read the whole graph again (a GFA of 5-115 MB with its
sequences) and search it again, about 60,000 times for the five graphs.
Simple-Read (our static, run the same way) does the same. Per batch they run
once per haplotype, a few hundred times, so it takes minutes.
Our algorithms take seconds to minutes in both modes.

BubbleFinder: the binary is expected at ~/tools/bubblefinder-1.0.5/BubbleFinder
(or set BUBBLEFINDER=/path/to/BubbleFinder); without it BubbleFinder is skipped
with a warning. Install (Linux x86-64):

    mkdir -p ~/tools/bubblefinder-1.0.5 && cd ~/tools/bubblefinder-1.0.5 \\
      && wget https://github.com/algbio/BubbleFinder/releases/download/v1.0.5/BubbleFinder-linux-x86_64.tar.gz \\
      && tar xzf BubbleFinder-linux-x86_64.tar.gz && cd -

=====================================================================
WHAT IT DOES
=====================================================================

  1. data      unpacks testdata/small/chromosome_<c>/*.gz into
               input/Human/chromosome_<c>/Graphs/ (base graph, haplotypes,
               edge list; built with minigraph/minichain from HPRC assemblies
               on T2T-CHM13, see the README)
  2. ours      scripts/superbubbleRunner.py (compiles and runs):
               per edge:  static (Simple), incRecompute (ImpDFS), incNormalBST,
                          incBST, incMinIBST (ImpDS), incMinIBSTH1 (ImpSB)
               per batch: static (Simple), incRecomputeBatch (ImpDFS),
                          incBSTBatch (ImpDS), incBSTBatchH2 (ImpSB)
               then Preprocessing/parse_to_csv.py (CSV and peak memory)
  3. tools     scripts/superbubbleComparator.py -t bubblegun / bubblefinder,
               per batch (per edge with --tools-per-edge); and Simple-Read
               (superbubbleRunner.py -a staticRead: static, but after every
               update it reads the graph from the same GFA as the tools and
               builds everything from scratch), in the same modes
  4. tables    output/artifactTest/small_tables.md (printed at the end):
               total time per graph and speedup over Simple, time of every
               step, minI / maxO of the data structures, the tools against
               our algorithms (time with the tools' reading and search
               apart, and memory), and whether everybody finds the same
               superbubbles in the final graph

Times are the sums of the timed steps of every update (our algorithms, without
reading the input) and reading the GFA plus the search (tools, after every
update), as in the paper. The tools' time is also split into reading and the
search alone; "search only" compares the same work on both sides. Times depend
on the machine, the speedups much less.
"""
import argparse
import os
import pathlib
import gzip
import shutil
import subprocess
import sys
import time

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "Preprocessing"))
from componentTable import steps, peak_kb, tool_peak_kb                                      # noqa: E402
from toolComparisonToLatex import (comparable, read_algorithm_superbubbles, read_final_graph,  # noqa: E402
                                   read_tool_superbubbles, tool_time_split, algorithm_time_split)
from superbubbleTools.bubblefinderTool import DEFAULT_BINARY                                  # noqa: E402

CHROMS = ["11", "17", "18", "22", "X"]
FILES = ["numeric_edgelist.txt", "reduced_graph.gfa", "ordered_haplotypes.txt", "unique_elements.txt"]
# (paper name, algorithm ID) per mode, as in the paper's tables
ALGOS = {"perEdge": [("Simple", "static"), ("ImpDFS", "incRecompute"), ("ImpDS", "incMinIBST"),
                     ("ImpSB", "incMinIBSTH1")],
         "perBatch": [("Simple", "static"), ("ImpDFS", "incRecomputeBatch"), ("ImpDS", "incBSTBatch"),
                      ("ImpSB", "incBSTBatchH2")]}
# per edge only: the data structures for minI / maxO (Simple, normal BST, BST, minI BST)
TREES = [("Simple", "static"), ("NBST", "incNormalBST"), ("BST", "incBST"), ("MinI", "incMinIBST")]
MODE_DIR = {"perEdge": "perEdge", "perBatch": "perbatch"}
TOOLS = [("BubbleGun", "bubblegun"), ("BubbleFinder", "bubblefinder")]
STEPS = ["Insert", "DFS", "Top. order", "MinI", "MaxO", "SB compute", "Total"]
OUT = REPO / "output" / "artifactTest"


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def run(cmd):
    """runs a step of the pipeline, its output goes to the terminal; stops on an error"""
    log("$ " + " ".join(cmd))
    if subprocess.run(cmd, cwd=REPO).returncode != 0:
        sys.exit(f"[!] failed: {' '.join(cmd)}")


def unpack(chroms):
    for c in chroms:
        src = REPO / "testdata" / "small" / f"chromosome_{c}"
        dst = REPO / "input" / "Human" / f"chromosome_{c}" / "Graphs"
        dst.mkdir(parents=True, exist_ok=True)
        for f in FILES:
            with gzip.open(src / f"{f}.gz", "rb") as fin, open(dst / f, "wb") as fout:
                shutil.copyfileobj(fin, fout)
        log(f"chr{c}: test data unpacked into {dst.relative_to(REPO)}")


def ds(c):
    return f"Human/chromosome_{c}"


def tool_split(c, tool, mode):
    """{"total", "read", "find"} (s) of a tool run: the whole time, reading the GFA, the search alone;
    None if it did not run"""
    return tool_time_split(REPO / "output" / ds(c) / tool / MODE_DIR[mode] / f"{tool}_time.txt")


def static_read_peak_kb(c, mode):
    """peak memory (KB) of staticRead, from the end of its time file"""
    for line in open(REPO / "output" / ds(c) / "timeFile" / MODE_DIR[mode] / "staticRead_time.txt"):
        if line.startswith("PeakMemory(KB):"):
            return float(line.split()[1])
    return float("nan")


def table(head, rows):
    """a Markdown table: head = column names, rows = lists of cells"""
    out = ["| " + " | ".join(head) + " |", "|" + "|".join(["---"] + ["---:"] * (len(head) - 1)) + "|"]
    out += ["| " + " | ".join(str(x) for x in r) + " |" for r in rows]
    return "\n".join(out)


def f2(x):
    return "-" if x is None else f"{x:.2f}"


def time_tables(chroms, mode):
    """total time per graph with the speedup over Simple, and every step summed over the graphs"""
    names = [f"chr{c}" for c in chroms] + ["Total"]
    tot, comp = {}, {}
    for name, a in ALGOS[mode]:
        s = [steps(REPO / "output" / ds(c) / "timeFile" / MODE_DIR[mode] / f"{a}_time.txt") for c in chroms]
        tot[name] = [x[6] for x in s] + [sum(x[6] for x in s)]
        comp[name] = [sum(x[i] for x in s) for i in range(7)]
    rows = []
    for name, _ in ALGOS[mode]:
        rows.append([name] + [f"{v:.3f}" for v in tot[name]])
        if name != "Simple":
            rows.append(["  speedup"] + [f2(b / v) for b, v in zip(tot["Simple"], tot[name])])
    t1 = table(["Algorithm"] + names, rows)
    rows = [[name] + [f"{v:.3f}" for v in comp[name]] + ([f2(comp["Simple"][6] / comp[name][6])] if name != "Simple"
                                                          else ["1.00"]) for name, _ in ALGOS[mode]]
    t2 = table(["Algorithm"] + STEPS + ["Speedup"], rows)
    return t1, t2


def trees_table(chroms):
    """per edge: minI and maxO of the data structures, summed over the graphs"""
    rows = []
    for name, a in TREES:
        s = [steps(REPO / "output" / ds(c) / "timeFile" / "perEdge" / f"{a}_time.txt") for c in chroms]
        rows.append([name, f"{sum(x[3] for x in s):.3f}", f"{sum(x[4] for x in s):.3f}"])
    return table(["Data structure", "MinI (s)", "MaxO (s)"], rows)


def tools_tables(chroms, mode, tools):
    """time and memory of the tools and our algorithms; speedup and ratio against BubbleGun.
    A tool's time is split into reading the GFA and the search alone; our times have no reading
    (the timed steps), so "search only" compares the same work on both sides"""
    names = [f"chr{c}" for c in chroms] + ["Total"]
    times, mems, split = {}, {}, {}
    for name, t in tools:
        v = [tool_split(c, t, mode) for c in chroms]
        if None in v:
            continue
        split[name] = {k: [x[k] for x in v] + [sum(x[k] for x in v)] for k in ("total", "read", "find")}
        times[name] = split[name]["total"]
        mems[name] = [tool_peak_kb(ds(c), MODE_DIR[mode], t) / 1024 for c in chroms]
    # Simple-Read: static reading the graph from the tools' GFA after every update, split like a tool
    v = [algorithm_time_split(REPO / "output" / ds(c) / "timeFile" / MODE_DIR[mode] / "staticRead_time.txt")
         for c in chroms]
    if None not in v:
        split["Simple-Read"] = {k: [x[k] for x in v] + [sum(x[k] for x in v)] for k in ("total", "read", "find")}
        times["Simple-Read"] = split["Simple-Read"]["total"]
        mems["Simple-Read"] = [static_read_peak_kb(c, mode) / 1024 for c in chroms]
    for name, a in ALGOS[mode]:
        v = [steps(REPO / "output" / ds(c) / "timeFile" / MODE_DIR[mode] / f"{a}_time.txt")[6] for c in chroms]
        times[name] = v + [sum(v)]
        mems[name] = [peak_kb(ds(c), MODE_DIR[mode], a)[0] / 1024 for c in chroms]
    base = "BubbleGun" if "BubbleGun" in times else next(iter(times))
    rows = []
    for name in times:
        rows.append([name] + [f"{v:.3f}" for v in times[name]])
        if name in split:
            rows.append(["  reading the GFA"] + [f"{r:.3f} ({100 * r / t:.0f} %)" for r, t in
                                                 zip(split[name]["read"], split[name]["total"])])
            rows.append(["  search only"] + [f"{v:.3f}" for v in split[name]["find"]])
        if name != base:
            rows.append([f"  speedup over {base}"] + [f2(b / v) for b, v in zip(times[base], times[name])])
            if base in split:
                ref = split[base]["find"]
                mine = split[name]["find"] if name in split else times[name]
                rows.append([f"  speedup over {base}, search only"] + [f2(b / v) for b, v in zip(ref, mine)])
    t1 = table(["Program"] + names, rows)
    rows = []
    for name in mems:
        rows.append([name] + [f"{v:.1f}" for v in mems[name]])
        if name != base:
            rows.append([f"  {base} / {name}"] + [f2(b / v) for b, v in zip(mems[base], mems[name])])
    t2 = table(["Program"] + names[:-1], rows)
    return t1, t2


def agreement_table(chroms, modes, tools):
    """final graph: do all our algorithms give the same superbubbles, and do the tools find them too?"""
    rows = []
    for c in chroms:
        graph = read_final_graph(REPO / "input" / ds(c) / "Graphs" / "numeric_edgelist.txt")
        lists = {}
        for mode in modes:
            for name, a in ALGOS[mode] + (TREES if mode == "perEdge" else []) + [("Simple-Read", "staticRead")]:
                p = read_algorithm_superbubbles(REPO / "output" / ds(c) / "timeFile" / MODE_DIR[mode] /
                                                f"{a}_superbubbles.txt")
                if p is not None:
                    lists[(mode, a)] = comparable(p, graph)
        ours = next(iter(lists.values()))
        row = [f"chr{c}", len(ours), "yes" if all(v == ours for v in lists.values()) else "NO"]
        for name, t in tools:
            theirs = None
            for mode in modes:
                p = read_tool_superbubbles(REPO / "output" / ds(c) / t / MODE_DIR[mode] / f"{t}_superbubbles.txt")
                if p is not None:
                    theirs = comparable(p, graph)
                    break
            row += ["-", "-", "-"] if theirs is None else [len(ours & theirs), len(ours - theirs), len(theirs - ours)]
        rows.append(row)
    head = ["Graph", "Superbubbles (ours)", "All our algorithms agree"]
    for name, _ in tools:
        head += [f"Both (ours, {name})", "Only ours", f"Only {name}"]
    return table(head, rows)


def main():
    ap = argparse.ArgumentParser(description="Artifact test: the small human graphs, all algorithms and tools")
    ap.add_argument("-c", "--chromosomes", nargs="+", choices=CHROMS, default=CHROMS)
    ap.add_argument("--runs", type=int, default=1)
    ap.add_argument("--tools-per-edge", action="store_true", help="also the tools and Simple-Read per edge (takes hours)")
    ap.add_argument("--no-tools", action="store_true", help="only our algorithms")
    ap.add_argument("--overwrite", action="store_true", help="replace results that exist in output/")
    args = ap.parse_args()
    chroms = [c for c in CHROMS if c in args.chromosomes]
    start = time.time()

    existing = [c for c in chroms if (REPO / "output" / ds(c) / "timeFile").exists()]
    if existing and not args.overwrite:
        sys.exit(f"[!] output/Human/chromosome_{{{','.join(existing)}}}/ has results already; they would be "
                 f"replaced. Use --overwrite to go on.")

    tools = [] if args.no_tools else list(TOOLS)
    bf = os.environ.get("BUBBLEFINDER", DEFAULT_BINARY)
    if ("BubbleFinder", "bubblefinder") in tools and not os.access(bf, os.X_OK):
        log(f"[!] BubbleFinder not found at {bf}: skipped (install: see the top of this file or the README)")
        tools.remove(("BubbleFinder", "bubblefinder"))
    tool_modes = ["perBatch"] + (["perEdge"] if args.tools_per_edge else [])
    if tools and args.tools_per_edge:
        log("[!] --tools-per-edge: the tools and Simple-Read run after every edge, this takes HOURS "
            "(about 8 h for all five)")

    log("1. test data")
    unpack(chroms)

    log("2. our algorithms")
    cnames = [f"chromosome_{c}" for c in chroms]
    py = sys.executable
    for mode in ("perEdge", "perBatch"):
        algos = [a for _, a in ALGOS[mode]] + ([a for _, a in TREES[1:3]] if mode == "perEdge" else [])
        run([py, "scripts/superbubbleRunner.py", "-d", "Human", "-c", *cnames, "-m", mode, "--runs", str(args.runs),
             "-a", *algos])
        run([py, "Preprocessing/parse_to_csv.py", "-d", "human", "-c", *chroms, "-m", mode.lower()])

    log("3. tools, and Simple-Read (static reading the tools' GFA after every update)")
    for _, t in tools:
        for mode in tool_modes:
            run([py, "scripts/superbubbleComparator.py", "-d", "Human", "-c", *chroms, "-m", mode, "-t", t])
    if tools:
        for mode in tool_modes:
            run([py, "scripts/superbubbleRunner.py", "-d", "Human", "-c", *cnames, "-m", mode, "--runs", "1",
                 "-a", "staticRead"])

    log("4. tables")
    OUT.mkdir(parents=True, exist_ok=True)
    parts = [f"# Small human graphs: {', '.join('chr' + c for c in chroms)}",
             f"Our algorithms: {args.runs} run(s) each, averaged. Times in seconds (sum of the timed steps of every "
             "update, without reading the input; tools: reading the GFA + search after every update, also shown "
             "apart: \"search only\" is the same work as our algorithms do), memory in MB (peak of the process).",
             "Simple = static, ImpDFS = incremental DFS, ImpDS = + data structures for minI/maxO, "
             "ImpSB = + superbubble update."]
    for mode, title in (("perEdge", "Per edge"), ("perBatch", "Per batch")):
        t1, t2 = time_tables(chroms, mode)
        parts += [f"## {title}: total time per graph", t1, f"## {title}: time of every step (all graphs)", t2]
        if mode == "perEdge":
            parts += ["## Per edge: minI and maxO with each data structure (all graphs)", trees_table(chroms)]
        if tools and mode in tool_modes:
            t1, t2 = tools_tables(chroms, mode, tools)
            parts += [f"## {title}: tools and our algorithms, time", t1,
                      f"## {title}: tools and our algorithms, peak memory (MB)", t2]
        elif tools:
            parts += [f"## {title}: tools", "Not run (python3 scripts/smallGraphsTest.py --tools-per-edge: takes hours)."]
    parts += ["## Superbubbles of the final graph (single-edge superbubbles and the artificial source/sink not "
              "counted)", agreement_table(chroms, ["perBatch", "perEdge"], tools)]
    text = "\n\n".join(parts) + "\n"
    (OUT / "small_tables.md").write_text(text)
    print("\n" + text)
    log(f"done in {(time.time() - start) / 60:.1f} min; tables in {(OUT / 'small_tables.md').relative_to(REPO)}")


if __name__ == "__main__":
    main()
