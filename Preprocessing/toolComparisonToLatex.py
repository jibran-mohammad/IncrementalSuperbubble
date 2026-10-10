"""
toolComparisonToLatex.py

Compares the best algorithm of this repository (BEST_ALGORITHM in
scripts/pipelineRunner.py, currently incMinIBST) with an external
superbubble tool (BubbleGun, BubbleFinder, ...) and writes LaTeX tables:
one with time / speedup / memory, one with how many superbubbles of the
final graph both found.

=====================================================================
HOW TO RUN
=====================================================================

Run all commands from the repository root, with the virtual environment
active (source venv/bin/activate).

Step 1 - produce the results this script reads (for every dataset and
mode you want in the tables; replace Ecoli/perEdge as needed):

    # the tool: time per update, peak memory, bubble lists
    #   -> output/Ecoli/bubblegun/perEdge/bubblegun_time.txt
    #      output/Ecoli/bubblegun/perEdge/bubblegun_superbubbles.txt
    python3 scripts/superbubbleComparator.py -d Ecoli -m perEdge -t bubblegun

    # the algorithm, end-to-end time (mean of 10 runs)
    #   -> output/Ecoli/incsuperbubble/perEdge/incMinIBST_time.txt
    python3 scripts/pipelineRunner.py -d Ecoli -m perEdge

    # the algorithm, peak memory and superbubble list of the final graph
    #   -> output/Ecoli/timeFile/perEdge/incMinIBST_time.txt
    #      output/Ecoli/timeFile/perEdge/incMinIBST_superbubbles.txt
    python3 scripts/superbubbleRunner.py -d Ecoli -m perEdge

    For a human chromosome use -d Human -c 14 (superbubbleComparator.py,
    pipelineRunner.py: -c c14) or -c chromosome_14 (superbubbleRunner.py).
    Results produced on a server can be copied into the same output/ folders.

Step 2 - make the tables:

    # E. coli, BubbleGun, both modes
    python3 Preprocessing/toolComparisonToLatex.py -d Ecoli

    # E. coli and chromosomes 11, 13 and X, only perEdge
    python3 Preprocessing/toolComparisonToLatex.py -d Ecoli Human -c 11 13 X -m perEdge

    # all chromosomes that have results of the tool
    python3 Preprocessing/toolComparisonToLatex.py -d Human

    # another tool (needs its results from superbubbleComparator.py -t bubblefinder)
    python3 Preprocessing/toolComparisonToLatex.py -d Ecoli -t bubblefinder

    # another algorithm instead of BEST_ALGORITHM
    python3 Preprocessing/toolComparisonToLatex.py -d Ecoli -a incBST

Options:

    -d / --datasets      Ecoli and/or Human                    (default: Ecoli)
    -c / --chromosomes   Human chromosomes, e.g. 11 13 X       (default: all
                         chromosomes that have a folder for the tool)
    -m / --modes         perEdge and/or perBatch               (default: both)
    -t / --tool          tool name = name of its output folder (default: bubblegun)
    --tool-label         tool name shown in the tables         (default: BubbleGun,
                         BubbleFinder, or the -t value)
    -a / --algorithm     algorithm ID, see ALGORITHMS.md       (default: BEST_ALGORITHM)
    --outdir             where the tables go                   (default:
                         output/toolComparison/<tool>/)

Missing results do not stop the script: the value is shown as "-" and a
line starting with "[!]" says which file is missing.

=====================================================================
OUTPUT (in output/toolComparison/<tool>/)
=====================================================================

    performance_perEdge.tex / performance_perBatch.tex   (+ .csv)
        One row per dataset:
        Updates                  number of updates (edges or haplotypes)
        <algorithm> time (s)     AlgorithmTime of pipelineRunner.py: reading
                                 the edge list + base graph + all updates,
                                 mean of its --repeat runs
        <tool> time (s)          TotalTime of the tool: for every update,
                                 reading the whole GFA + finding bubbles
        Speedup                  tool time / algorithm time
        <algorithm> memory (MB)  PeakMemory(KB) of superbubbleRunner.py
        <tool> memory (MB)       largest PeakMemory(KB) of the tool
        Preprocessing (s)        GFA -> input of both (edge list, reduced
                                 GFA); shared by both, so not added to either

    agreement.tex   (+ .csv)
        One row per dataset, final graph: number of superbubbles of the
        algorithm, of the tool, found by both, only by the algorithm, only
        by the tool. If the algorithm is correct and the tool reports all
        superbubbles, both "Only" columns are 0.

    agreement_details_<dataset>.txt
        The superbubbles counted in the two "Only" columns, by GFA segment
        name, to look up in the GFA.

How superbubbles are compared: as unordered pairs of GFA segment names
(input id k = segment s{k-1}, the edge-list convention). Left out on both
sides: pairs with the artificial source/sink (ids 0 and 1), which only
the algorithm has, and single-edge superbubbles (the entrance's only
out-neighbour is the exit), which tools like BubbleGun do not report.
Both sides are filtered with the same final graph, read from
input/<dataset>/Graphs/numeric_edgelist.txt.

=====================================================================
ADDING ANOTHER TOOL (e.g. BubbleFinder)
=====================================================================

Nothing changes in this script. Write a tool class in
Preprocessing/superbubbleTools/ (name = "bubblefinder", run(), memory();
see toolInterface.py and bubblegunTool.py), register it in TOOLS in
scripts/superbubbleComparator.py, run that with -t bubblefinder, then run
this script with -t bubblefinder. Add the tool to TOOL_LABELS below if
the name in the tables should be spelled differently.
"""

import argparse
import collections
import pathlib
import sys

import pandas as pd

# ------------------------------------------------------------------
# Setup: find the repository and reuse settings from other scripts
# ------------------------------------------------------------------

# Repository root = parent of this file's directory (Preprocessing/..), so
# the script works whatever the repository folder is called.
REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent

# Make "Preprocessing.csvToLatex" and "pipelineRunner" importable.
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "scripts"))

# DISPLAY_NAMES: algorithm ID -> name shown in tables (e.g. Inc-MinI-BST);
# format_float / generate_latex_table: same number format and table layout
# as the other LaTeX tables of this repository.
from Preprocessing.csvToLatex import DISPLAY_NAMES, format_float, generate_latex_table

# BEST_ALGORITHM: default algorithm; MODE_DIRS: perEdge -> perEdge,
# perBatch -> perbatch (folder names); PIPELINE_DIR: folder of the
# pipelineRunner results (incsuperbubble).
from pipelineRunner import BEST_ALGORITHM, MODE_DIRS, TOOL_NAME as PIPELINE_DIR


# Names used in table headers; any other tool is shown as given with -t
# (or set --tool-label).
TOOL_LABELS = {
    "bubblegun": "BubbleGun",
    "bubblefinder": "BubbleFinder",
}


############################################################
# Reading result files
#
# Every reader returns an empty result ({} or None) when the file does not
# exist, so a missing file shows up as "-" in the table instead of an error.
############################################################

def read_pipeline_time(path):
    """AlgorithmTime / AlgorithmTimeStd / PreprocessingTime / Runs from pipelineRunner."""

    values = {}

    if not path.exists():
        return values

    # Lines look like "AlgorithmTime: 0.261753 seconds" or "Runs: 10".
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

    # Last line of the timing file: "PeakMemory(KB): 5105" (may be a
    # decimal number, because superbubbleRunner averages 100 runs).
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

    # CSV written by superbubbleComparator.py: one row per update whose
    # first field is the step number, a "Base" row, and summary lines
    # "TotalTime,<s>" and "PeakMemory(KB),<KB>" at the end.
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

    # Superbubble files consist of blocks, each starting with a "# ..."
    # header (tool: base graph and after every haplotype; algorithm: only
    # the final graph). Every header starts a new block, so after the loop
    # `block` holds the last one = the final graph, split into fields.
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

    # e.g. "simple s3157 s3155": keep the two ends, drop the bubble type.
    return [(fields[1], fields[2]) for fields in last_block(path)]


def read_algorithm_superbubbles(path):
    """Final-graph superbubbles of the algorithm: lines '<entrance> <exit>' (input ids)."""

    if not path.exists():
        return None

    # e.g. "927 928" -> ("s926", "s927"): input id k is GFA segment s{k-1}.
    # 0 and 1 (artificial source and sink) have no segment and are kept as
    # "0"/"1", so comparable() can drop them.
    return [(f"s{int(a) - 1}" if int(a) > 1 else a, f"s{int(b) - 1}" if int(b) > 1 else b)
            for a, b in (fields[:2] for fields in last_block(path))]


def read_final_graph(path):
    """Out-neighbours (GFA names) of the final graph: base edges + all updates."""

    out = collections.defaultdict(set)

    if not path.exists():
        return None

    # numeric_edgelist.txt: comment lines, the graph name, a header line
    # "<number of base edges> <number of nodes>", the base edges "u v",
    # then for every haplotype a "# HAPLOTYPE: ..." line, its number of
    # edges and its edges "u v". Every line with exactly two numbers is an
    # edge, except the first one, which is the header and must not be read
    # as an edge (it would add a fake edge between two node ids).
    header_seen = False

    for line in open(path):
        fields = line.split()
        if len(fields) == 2 and fields[0].isdigit() and fields[1].isdigit():
            if not header_seen:
                header_seen = True
                continue
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

        # The artificial source (0) and sink (1) exist only in the algorithm.
        if a in ("0", "1") or b in ("0", "1"):
            continue

        # Single-edge superbubble: nothing in between entrance and exit.
        # Checked in both directions, because a tool may report the ends
        # in either order.
        if out.get(a) == {b} or out.get(b) == {a}:
            continue

        # frozenset = unordered pair, so (s1, s2) and (s2, s1) are equal.
        result.add(frozenset((a, b)))

    return result


############################################################
# Building the table rows
############################################################

def dataset_label(dataset_path):
    """Ecoli -> Ecoli, Human/chromosome_11 -> chr11"""

    name = dataset_path.name
    return "chr" + name.split("_", 1)[1] if name.startswith("chromosome_") else name


def performance_rows(datasets, mode, algorithm, tool):
    """One row per dataset for performance_<mode>: times, speedup, memory."""

    rows = []

    for ds in datasets:

        out = REPO_ROOT / "output" / ds
        mode_dir = MODE_DIRS[mode]

        # Algorithm: end-to-end time (pipelineRunner) and peak memory
        # (superbubbleRunner); tool: time, updates and peak memory.
        ours = read_pipeline_time(out / PIPELINE_DIR / mode_dir / f"{algorithm}_time.txt")
        ours_memory = read_algorithm_memory(out / "timeFile" / mode_dir / f"{algorithm}_time.txt")
        theirs = read_tool_time(out / tool / mode_dir / f"{tool}_time.txt")

        # Dataset not run in this mode at all: no row.
        if not ours and not theirs:
            continue

        # Partly missing: still a row (with "-"), plus a message saying
        # which script has to be run.
        for what, ok in (("algorithm time (pipelineRunner.py)", ours),
                         ("algorithm memory (superbubbleRunner.py)", ours_memory),
                         (f"{tool} results (superbubbleComparator.py)", theirs)):
            if not ok:
                print(f"  [!] {ds} {mode}: no {what}")

        t_ours = ours.get("AlgorithmTime")
        t_theirs = theirs.get("TotalTime")

        # Times in seconds, memory converted from KB to MB; "-" if missing.
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
    """One row per dataset for agreement: superbubbles of the final graph in both / only one."""

    rows = []

    for ds in datasets:

        out = REPO_ROOT / "output" / ds

        # Final graph, used to recognise single-edge superbubbles.
        graph = read_final_graph(REPO_ROOT / "input" / ds / "Graphs" / "numeric_edgelist.txt")

        # The final graph is the same in both modes; use the first mode
        # for which both lists exist.
        for mode in modes:
            mode_dir = MODE_DIRS[mode]
            ours = read_algorithm_superbubbles(out / "timeFile" / mode_dir / f"{algorithm}_superbubbles.txt")
            theirs = read_tool_superbubbles(out / tool / mode_dir / f"{tool}_superbubbles.txt")
            if ours is not None and theirs is not None:
                break

        # Something missing: no row. Say what is missing, unless the dataset
        # has no lists at all (then it simply was not run).
        if graph is None or ours is None or theirs is None:
            if ours is not None or theirs is not None:
                missing = [what for what, v in (("final graph (numeric_edgelist.txt)", graph),
                                                (f"{algorithm}_superbubbles.txt", ours),
                                                (f"{tool}_superbubbles.txt", theirs)) if v is None]
                print(f"  [!] {ds}: no agreement row, missing {', '.join(missing)}")
            continue

        # Same filtering on both sides, then compare as sets.
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

        # The superbubbles behind the two "Only" columns, for inspection.
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

        # Ecoli (or any other non-human dataset): one folder output/<name>/.
        if name != "Human":
            datasets.append(pathlib.Path(name))
            continue

        # Human: the chromosomes given with -c ("14" or "chromosome_14"),
        # otherwise every chromosome that has a results folder for the tool.
        if chromosomes:
            chrs = [c if c.startswith("chromosome_") else f"chromosome_{c}" for c in chromosomes]
        else:
            human = REPO_ROOT / "output" / "Human"
            chrs = sorted(d.name for d in human.glob("chromosome_*") if (d / tool).is_dir())

        datasets += [pathlib.Path("Human") / c for c in chrs]

    return datasets


############################################################
# Main: arguments, then one performance table per mode and
# one agreement table
############################################################

def main():

    parser = argparse.ArgumentParser(
        description="Compare the best algorithm with a superbubble tool (BubbleGun, "
                    "BubbleFinder, ...) and write LaTeX tables. See the top of this file "
                    "for which results to produce first.")

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

    # Names shown in the tables, e.g. "Inc-MinI-BST" and "BubbleGun".
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

        # Internal row keys -> column headers of the table.
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

        # Same table as CSV (easy to check or to use elsewhere) and as LaTeX.
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

        # Internal row keys -> column headers ("Mode" stays in the CSV only).
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
