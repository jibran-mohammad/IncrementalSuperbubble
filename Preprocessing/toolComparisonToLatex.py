"""
toolComparisonToLatex.py

Compares the best algorithm of this repository (BEST_ALGORITHM in
scripts/pipelineRunner.py, e.g. incMinIBST) with an external superbubble
tool (BubbleGun, BubbleFinder, ...) and writes LaTeX tables.

The tool is chosen with -t and only needs the two files that
scripts/superbubbleComparator.py writes for every tool:

    output/<dataset>/<tool>/<perEdge|perbatch>/<tool>_time.txt
    output/<dataset>/<tool>/<perEdge|perbatch>/<tool>_superbubbles.txt

so a new tool needs no change here: add it to superbubbleComparator.py,
run that, then run this script with -t <tool>.

Inputs for the algorithm:

    output/<dataset>/incsuperbubble/<mode>/<ID>_time.txt    pipelineRunner.py:
        AlgorithmTime (mean of --repeat runs: reading the edge list + base
        graph + all updates), PreprocessingTime
    output/<dataset>/timeFile/<mode>/<ID>_time.txt          superbubbleRunner.py:
        PeakMemory(KB)
    output/<dataset>/timeFile/<mode>/<ID>_superbubbles.txt  superbubbleRunner.py:
        superbubbles of the final graph (input ids)

Tables (in --outdir, default output/toolComparison/<tool>/):

    performance_<mode>.tex   per dataset: number of updates, algorithm time,
                             tool time, speedup, peak memory of both,
                             preprocessing time (shared by both: the tool's
                             input comes from the same preprocessing, so it
                             is reported separately and not added to either)
    agreement.tex            per dataset, final graph: superbubbles found by
                             both / only the algorithm / only the tool
    agreement_details_<dataset>.txt   the differing superbubbles

Superbubbles are compared as unordered pairs of GFA segment names (input
id k = segment s{k-1}). Pairs with the artificial source/sink (0, 1) and
single-edge superbubbles (entrance whose only out-neighbour is the exit;
tools such as BubbleGun do not report them) are left out on both sides,
using the final graph from input/<dataset>/Graphs/numeric_edgelist.txt.

Usage:

    python3 Preprocessing/toolComparisonToLatex.py -d Ecoli
    python3 Preprocessing/toolComparisonToLatex.py -d Ecoli Human -c 11 13 X -t bubblegun
    python3 Preprocessing/toolComparisonToLatex.py -d Human -t bubblefinder -m perEdge
"""

import argparse
import collections
import pathlib
import sys

import pandas as pd

# Repository root = parent of this file's directory (no dependency on the
# name of the repository folder).
REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent

sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from Preprocessing.csvToLatex import DISPLAY_NAMES, format_float, generate_latex_table
from pipelineRunner import BEST_ALGORITHM, MODE_DIRS, TOOL_NAME as PIPELINE_DIR


# Names used in table headers; any other tool is shown as given with -t
# (or set --tool-label).
TOOL_LABELS = {
    "bubblegun": "BubbleGun",
    "bubblefinder": "BubbleFinder",
}


############################################################
# Reading result files
############################################################

def read_pipeline_time(path):
    """AlgorithmTime / AlgorithmTimeStd / PreprocessingTime / Runs from pipelineRunner."""

    values = {}

    if not path.exists():
        return values

    for line in open(path):
        key, _, rest = line.partition(":")
        if key in ("AlgorithmTime", "AlgorithmTimeStd", "PreprocessingTime"):
            values[key] = float(rest.split()[0])
        elif key == "Runs":
            values[key] = int(rest)

    return values


def read_algorithm_memory(path):
    """PeakMemory(KB) from a superbubbleRunner timing file."""

    if not path.exists():
        return None

    for line in open(path):
        if line.startswith("PeakMemory(KB):"):
            return int(float(line.split(":")[1]))

    return None


def read_tool_time(path):
    """Number of updates, TotalTime and PeakMemory(KB) from <tool>_time.txt."""

    values = {}

    if not path.exists():
        return values

    updates = 0

    for line in open(path):
        fields = line.strip().split(",")
        if fields[0].isdigit():
            updates += 1
        elif fields[0] == "TotalTime":
            values["TotalTime"] = float(fields[1])
        elif fields[0] == "PeakMemory(KB)":
            values["PeakMemory"] = int(fields[1])

    values["Updates"] = updates

    return values


def last_block(path):
    """Lines of the last '#'-headed block of a superbubble file (the final graph)."""

    block = None

    for line in open(path):
        line = line.strip()
        if line.startswith("#"):
            block = []
        elif line and block is not None:
            block.append(line.split())

    return block or []


def read_tool_superbubbles(path):
    """Final-graph bubbles of the tool: lines '<type> <end1> <end2>' (GFA names)."""

    if not path.exists():
        return None

    return [(fields[1], fields[2]) for fields in last_block(path)]


def read_algorithm_superbubbles(path):
    """Final-graph superbubbles of the algorithm: lines '<entrance> <exit>' (input ids)."""

    if not path.exists():
        return None

    return [(f"s{int(a) - 1}" if int(a) > 1 else a, f"s{int(b) - 1}" if int(b) > 1 else b)
            for a, b in (fields[:2] for fields in last_block(path))]


def read_final_graph(path):
    """Out-neighbours (GFA names) of the final graph: base edges + all updates."""

    out = collections.defaultdict(set)

    if not path.exists():
        return None

    for line in open(path):
        fields = line.split()
        if len(fields) == 2 and fields[0].isdigit() and fields[1].isdigit():
            out[f"s{int(fields[0]) - 1}"].add(f"s{int(fields[1]) - 1}")

    return out


############################################################
# Comparing superbubbles
############################################################

def comparable(pairs, out):
    """
    Unordered pairs without the artificial source/sink (ids 0, 1) and without
    single-edge superbubbles (entrance whose only out-neighbour is the exit).
    """

    result = set()

    for a, b in pairs:

        if a in ("0", "1") or b in ("0", "1"):
            continue

        if out.get(a) == {b} or out.get(b) == {a}:
            continue

        result.add(frozenset((a, b)))

    return result


############################################################
# Tables
############################################################

def dataset_label(dataset_path):
    """Ecoli -> Ecoli, Human/chromosome_11 -> chr11"""

    name = dataset_path.name
    return "chr" + name.split("_", 1)[1] if name.startswith("chromosome_") else name


def performance_rows(datasets, mode, algorithm, tool):

    rows = []

    for ds in datasets:

        out = REPO_ROOT / "output" / ds
        mode_dir = MODE_DIRS[mode]

        ours = read_pipeline_time(out / PIPELINE_DIR / mode_dir / f"{algorithm}_time.txt")
        ours_memory = read_algorithm_memory(out / "timeFile" / mode_dir / f"{algorithm}_time.txt")
        theirs = read_tool_time(out / tool / mode_dir / f"{tool}_time.txt")

        if not ours and not theirs:
            continue

        for what, ok in (("algorithm time (pipelineRunner.py)", ours),
                         ("algorithm memory (superbubbleRunner.py)", ours_memory),
                         (f"{tool} results (superbubbleComparator.py)", theirs)):
            if not ok:
                print(f"  [!] {ds} {mode}: no {what}")

        t_ours = ours.get("AlgorithmTime")
        t_theirs = theirs.get("TotalTime")

        rows.append({
            "Dataset": dataset_label(ds),
            "Updates": theirs.get("Updates", "-"),
            "AlgoTime": format_float(t_ours) if t_ours is not None else "-",
            "ToolTime": format_float(t_theirs) if t_theirs is not None else "-",
            "Speedup": format_float(t_theirs / t_ours, precision=1) if t_ours and t_theirs else "-",
            "AlgoMem": format_float(ours_memory / 1024, precision=1) if ours_memory else "-",
            "ToolMem": format_float(theirs["PeakMemory"] / 1024, precision=1) if "PeakMemory" in theirs else "-",
            "Preprocessing": format_float(ours["PreprocessingTime"]) if "PreprocessingTime" in ours else "-",
        })

    return rows


def agreement_rows(datasets, modes, algorithm, tool, outdir):

    rows = []

    for ds in datasets:

        out = REPO_ROOT / "output" / ds

        graph = read_final_graph(REPO_ROOT / "input" / ds / "Graphs" / "numeric_edgelist.txt")

        # The final graph is the same in both modes; use the first mode
        # for which both lists exist.
        for mode in modes:
            mode_dir = MODE_DIRS[mode]
            ours = read_algorithm_superbubbles(out / "timeFile" / mode_dir / f"{algorithm}_superbubbles.txt")
            theirs = read_tool_superbubbles(out / tool / mode_dir / f"{tool}_superbubbles.txt")
            if ours is not None and theirs is not None:
                break

        if graph is None or ours is None or theirs is None:
            if ours is not None or theirs is not None:
                missing = [what for what, v in (("final graph (numeric_edgelist.txt)", graph),
                                                (f"{algorithm}_superbubbles.txt", ours),
                                                (f"{tool}_superbubbles.txt", theirs)) if v is None]
                print(f"  [!] {ds}: no agreement row, missing {', '.join(missing)}")
            continue

        a = comparable(ours, graph)
        b = comparable(theirs, graph)

        rows.append({
            "Dataset": dataset_label(ds),
            "Mode": mode,
            "Algo": len(a),
            "Tool": len(b),
            "Both": len(a & b),
            "OnlyAlgo": len(a - b),
            "OnlyTool": len(b - a),
        })

        details = outdir / f"agreement_details_{dataset_label(ds)}.txt"

        with open(details, "w") as f:
            f.write(f"# {ds}, final graph ({mode}): {algorithm} vs {tool}\n")
            f.write(f"# only {algorithm}: {len(a - b)}\n")
            for p in sorted(sorted(p) for p in a - b):
                f.write(f"{p[0]} {p[1]}\n")
            f.write(f"# only {tool}: {len(b - a)}\n")
            for p in sorted(sorted(p) for p in b - a):
                f.write(f"{p[0]} {p[1]}\n")

    return rows


def discover_datasets(names, chromosomes, tool):
    """Ecoli -> [Ecoli]; Human -> Human/chromosome_X for -c or all with tool output."""

    datasets = []

    for name in names:

        if name != "Human":
            datasets.append(pathlib.Path(name))
            continue

        if chromosomes:
            chrs = [c if c.startswith("chromosome_") else f"chromosome_{c}" for c in chromosomes]
        else:
            human = REPO_ROOT / "output" / "Human"
            chrs = sorted(d.name for d in human.glob("chromosome_*") if (d / tool).is_dir())

        datasets += [pathlib.Path("Human") / c for c in chrs]

    return datasets


############################################################

def main():

    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[1])

    parser.add_argument("-d", "--datasets", nargs="+", default=["Ecoli"],
                        help="Ecoli and/or Human (default: Ecoli)")

    parser.add_argument("-c", "--chromosomes", nargs="+",
                        help="Human chromosomes, e.g. 11 13 X (default: all with tool output)")

    parser.add_argument("-m", "--modes", nargs="+", choices=["perEdge", "perBatch"],
                        default=["perEdge", "perBatch"])

    parser.add_argument("-t", "--tool", default="bubblegun",
                        help="Tool name = its output folder, e.g. bubblegun, bubblefinder")

    parser.add_argument("--tool-label",
                        help="Name of the tool in table headers (default: built-in label or -t)")

    parser.add_argument("-a", "--algorithm", default=BEST_ALGORITHM,
                        help=f"Algorithm ID (default: BEST_ALGORITHM = {BEST_ALGORITHM})")

    parser.add_argument("--outdir",
                        help="Output directory (default: output/toolComparison/<tool>)")

    args = parser.parse_args()

    tool_label = args.tool_label or TOOL_LABELS.get(args.tool, args.tool)
    algo_label = DISPLAY_NAMES.get(args.algorithm, args.algorithm)

    outdir = pathlib.Path(args.outdir) if args.outdir else REPO_ROOT / "output" / "toolComparison" / args.tool
    outdir.mkdir(parents=True, exist_ok=True)

    datasets = discover_datasets(args.datasets, args.chromosomes, args.tool)

    print(f"Comparing {args.algorithm} with {args.tool} on: {', '.join(map(str, datasets))}")

    ################ performance, one table per mode ################

    for mode in args.modes:

        rows = performance_rows(datasets, mode, args.algorithm, args.tool)

        if not rows:
            print(f"  [!] {mode}: no results, table skipped")
            continue

        df = pd.DataFrame(rows)

        columns = {
            "Dataset": "Dataset",
            "Updates": "Updates",
            "AlgoTime": f"{algo_label} time (s)",
            "ToolTime": f"{tool_label} time (s)",
            "Speedup": "Speedup",
            "AlgoMem": f"{algo_label} memory (MB)",
            "ToolMem": f"{tool_label} memory (MB)",
            "Preprocessing": "Preprocessing (s)",
        }

        df = df.rename(columns=columns)

        df.to_csv(outdir / f"performance_{mode}.csv", index=False)

        generate_latex_table(
            df, list(columns.values()),
            caption=(f"{algo_label} vs {tool_label} ({mode}): total time over all updates "
                     f"(both read their input and process every update), speedup = "
                     f"{tool_label} time / {algo_label} time, and peak memory. "
                     f"Preprocessing produces the input of both and is not included in either time."),
            label=f"tab:{args.tool}_performance_{mode}",
            out_tex=str(outdir / f"performance_{mode}.tex"),
        )

    ################ agreement on the final graph ################

    rows = agreement_rows(datasets, args.modes, args.algorithm, args.tool, outdir)

    if rows:

        df = pd.DataFrame(rows)

        columns = {
            "Dataset": "Dataset",
            "Algo": algo_label,
            "Tool": tool_label,
            "Both": "Both",
            "OnlyAlgo": f"Only {algo_label}",
            "OnlyTool": f"Only {tool_label}",
        }

        df = df.rename(columns=columns)

        df.to_csv(outdir / "agreement.csv", index=False)

        generate_latex_table(
            df, list(columns.values()),
            caption=(f"Superbubbles of the final graph found by {algo_label} and {tool_label} "
                     f"(unordered entrance/exit pairs; single-edge superbubbles and the "
                     f"artificial source/sink excluded)."),
            label=f"tab:{args.tool}_agreement",
            out_tex=str(outdir / "agreement.tex"),
        )

    else:
        print("  [!] no agreement table: superbubble lists missing "
              f"({args.algorithm}_superbubbles.txt / {args.tool}_superbubbles.txt)")

    print(f"Output: {outdir}")


if __name__ == "__main__":
    main()
