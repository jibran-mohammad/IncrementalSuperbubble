"""
superbubbleComparator.py

Runs an external superbubble tool (BubbleGun, later e.g. BubbleFinder) on
the same sequence of graph updates that the C++ algorithms process: after
every update (an edge in perEdge mode, a haplotype in perBatch mode) the
tool computes all bubbles of the whole graph from scratch. Records the
time of every run, the tool's peak memory and its bubble lists, for the
comparison with the algorithms (Preprocessing/toolComparisonToLatex.py).
What is timed and how the comparison works: COMPARISON_METHOD.md.

=====================================================================
HOW TO RUN
=====================================================================

Run from the repository root, with the virtual environment active
(source venv/bin/activate; BubbleGun is installed there). The graph files
must exist first (scripts/graphRunner.py).

    # BubbleGun on E. coli, one run per edge (about 1.5 minutes)
    python3 scripts/superbubbleComparator.py -d Ecoli -m perEdge

    # one run per haplotype
    python3 scripts/superbubbleComparator.py -d Ecoli -m perBatch

    # human chromosomes 14 and X
    python3 scripts/superbubbleComparator.py -d Human -c 14 X -m perEdge

    # every human chromosome that has been preprocessed
    python3 scripts/superbubbleComparator.py -d Human -m perEdge

    # E. coli and every human chromosome, one after the other
    python3 scripts/superbubbleComparator.py -d all -m perEdge

    # another tool, once it is added (see below)
    python3 scripts/superbubbleComparator.py -d Ecoli -m perEdge -t bubblefinder

Options:

    -d           Ecoli | Human | HPRC | all (= Ecoli + Human)    (required;
                 HPRC: Minigraph-Cactus graphs prepared by
                 Preprocessing/mcToDag.py, experimental)
    -c           Human only: chromosomes, e.g. 14 X or chromosome_14
                 (default: every chromosome with a reduced_graph.gfa)
    -m           perEdge | perBatch                              (required)
    -t / --tool  tool to run, by its name                        (default: bubblegun)

A dataset stops with an error if the tool reports a warning (for
BubbleGun e.g. a link to a segment that is not in the GFA), so wrong
input does not go unnoticed; the other datasets still run, and a SUMMARY
at the end lists every dataset as ok or FAILED (exit code 1 if any
failed).

=====================================================================
INPUT
=====================================================================

From input/Ecoli/Graphs/ or input/Human/chromosome_<X>/Graphs/:

    reduced_graph.gfa        base graph (reference only)
    ordered_haplotypes.txt   order in which haplotypes are added
    unique_elements.txt      nodes and edges each haplotype adds

The graph after each update (base graph + the edges added so far, with an
S line for every new node) is written to a temporary GFA and given to the
tool; nodes and edges are written in the order they were added, so the
temporary GFA, and the tool's output, is the same in every run.

=====================================================================
OUTPUT (folder named after the tool)
=====================================================================

    output/Ecoli/<tool>/<perEdge|perbatch>/
    output/Human/chromosome_<X>/<tool>/<perEdge|perbatch>/

    <tool>_time.txt          CSV, one row per update (edge or haplotype):
                             Time        one tool run on the whole graph
                                         (BubbleGun: reading the GFA +
                                         finding bubbles; Python start-up
                                         not included)
                             ReadTime / FindTime / ChainsTime   its parts
                                         (ChainsTime: building chains and
                                         nesting, not included in Time)
                             simple / super / insertion / TotalBubbles
                                         bubbles of each type (all three
                                         are superbubbles)
                             PeakMemoryKB / BaselineMemoryKB   peak memory
                                         of the tool, measured at the base
                                         graph and after each haplotype in
                                         a separate, untimed run ("-" on
                                         other rows); baseline = before the
                                         graph is read
                             then TotalUpdateTime, TotalTime, and the
                             largest PeakMemory(KB) / BaselineMemory(KB)
    <tool>_superbubbles.txt  every bubble ("<type> <end1> <end2>", GFA
                             segment names) for the base graph and after
                             each haplotype; the last block is the final
                             graph

=====================================================================
ADDING ANOTHER TOOL (e.g. BubbleFinder)
=====================================================================

  1. Write Preprocessing/superbubbleTools/bubblefinderTool.py with a class
     like BubbleGunTool (Preprocessing/superbubbleTools/bubblegunTool.py):
     name = "bubblefinder", run(gfa_file, keep_bubbles) and optionally
     memory(gfa_file); the expected return values are described in
     Preprocessing/superbubbleTools/toolInterface.py.
  2. Import it below and add it to TOOLS.
  3. Run this script with -t bubblefinder; results go to
     output/<dataset>/bubblefinder/.
"""

import sys
import pathlib
import argparse
import time


############################################################
# Dynamic Path Discovery
############################################################

def find_repo_root(start_path: pathlib.Path, repo_name: str) -> pathlib.Path:

    current_dir = start_path.resolve()

    while current_dir != current_dir.parent:

        if current_dir.name == repo_name:
            return current_dir

        current_dir = current_dir.parent

    if current_dir.name == repo_name:
        return current_dir

    raise FileNotFoundError(
        f"Could not find repository root '{repo_name}'"
    )


############################################################
# Repository Setup
############################################################

REPO_NAME = "IncSuperbubble"

try:

    try:
        script_dir = pathlib.Path(__file__).resolve().parent
    except NameError:
        script_dir = pathlib.Path.cwd()

    REPO_ROOT_PATH = find_repo_root(script_dir, REPO_NAME)

    print(f"ℹ️ Located repository root: {REPO_ROOT_PATH}")
except FileNotFoundError as e:

    print(e, file=sys.stderr)
    sys.exit(1)

############################################################
# Critical Line (Same Logic as graphRunner)
############################################################

sys.path.insert(0, str(REPO_ROOT_PATH))


############################################################
# Now imports work
############################################################


from Preprocessing.haplotypeParser import *
from Preprocessing.graphUpdater import *
from Preprocessing.superbubbleTools.bubblegunTool import BubbleGunTool


# Tools selectable with -t; the key is the tool's own name.
TOOLS = {tool.name: tool for tool in (BubbleGunTool,)}

# Same mode directory names as superbubbleRunner.py
MODE_DIRS = {"perEdge": "perEdge", "perBatch": "perbatch"}


#######################################################

def write_time_row(timefile,step,hap,update,result,memory=None):

    counts = result["counts"]

    phases = result.get("phases", {})

    # memory is only measured at the base graph and at haplotype boundaries
    mem = f"{memory['peak']},{memory['baseline']}" if memory else "-,-"

    timefile.write(
        f"{step},{hap},{update},{result['time']},"
        + "".join(f"{t}," for t in phases.values())
        + ",".join(str(n) for n in counts.values())
        + f",{sum(counts.values())},{mem}\n"
    )


def measure_memory(tool,gfa,memories):
    """Peak memory of the tool on gfa (separate, untimed run), if supported."""

    memory = tool.memory(gfa)

    if memory:
        memories.append(memory)

    return memory


def write_superbubbles(sbfile,title,result):

    counts = result["counts"]

    summary = " ".join(f"{t}={n}" for t,n in counts.items())

    sbfile.write(f"# {title} | {summary} total={sum(counts.values())}\n")

    for bubble_type,end1,end2 in sorted(result["bubbles"]):
        sbfile.write(f"{bubble_type} {end1} {end2}\n")

    sbfile.write("\n")


def run_updates(mode,updater,hap_order,hap_edges,tool,timefile,sbfile):

    ################################################
    # Base graph
    ################################################

    memories = []

    gfa = updater.write_temp_gfa()

    result = tool.run(gfa, keep_bubbles=True)

    base_time = result["time"]

    timefile.write(
        "Step,Haplotype," + ("Edge" if mode == "perEdge" else "Edges")
        + ",Time,"
        + "".join(f"{p.capitalize()}Time," for p in result.get("phases", {}))
        + ",".join(result["counts"]) + ",TotalBubbles,PeakMemoryKB,BaselineMemoryKB\n"
    )

    write_time_row(timefile,"Base","-","-",result,measure_memory(tool,gfa,memories))

    write_superbubbles(sbfile,"Base graph",result)

    print("Base completed")


    ################################################
    # Updates: one tool run per edge (perEdge) or
    # per haplotype (perBatch); bubbles are listed
    # after the last update of every haplotype
    ################################################

    update_time = 0

    step = 0

    for h, hap in enumerate(hap_order, 1):

        edges = hap_edges[hap]

        if mode == "perEdge":

            for e, edge in enumerate(edges, 1):

                updater.add_edge(edge)

                step += 1

                last = (e == len(edges))

                gfa = updater.write_temp_gfa()

                result = tool.run(gfa, keep_bubbles=last)

                update_time += result["time"]

                memory = measure_memory(tool,gfa,memories) if last else None

                write_time_row(timefile,step,hap,edge,result,memory)

        else:

            for edge in edges:
                updater.add_edge(edge)

            step += 1

            gfa = updater.write_temp_gfa()

            result = tool.run(gfa, keep_bubbles=True)

            update_time += result["time"]

            write_time_row(timefile,step,hap,len(edges),result,measure_memory(tool,gfa,memories))

        write_superbubbles(
            sbfile,
            f"After haplotype {h}/{len(hap_order)}: {hap} (step {step})",
            result
        )

        print(f"Haplotype {h}/{len(hap_order)} completed: {hap}")


    ################################################
    # Total tool runtime
    ################################################

    timefile.write("\n")

    timefile.write(f"TotalUpdateTime,{update_time}\n")

    timefile.write(f"TotalTime,{base_time + update_time}\n")

    print(f"Total {tool.name} time:", base_time + update_time)

    # Largest measured peak, comparable with PeakMemory(KB) of the C++ algorithms
    if memories:

        largest = max(memories, key=lambda m: m["peak"])

        timefile.write(f"PeakMemory(KB),{largest['peak']}\n")

        timefile.write(f"BaselineMemory(KB),{largest['baseline']}\n")

        print(f"{tool.name} peak memory: {largest['peak']} KB (baseline {largest['baseline']} KB)")


#######################################################

def run_dataset(repo, dataset_path, mode, tool_name):
    """Run the tool after every update of one dataset; write its two output files."""

    graph_dir = repo / "input" / dataset_path / "Graphs"

    base_gfa=graph_dir/"reduced_graph.gfa"

    order_file=graph_dir/"ordered_haplotypes.txt"

    unique_file=graph_dir/"unique_elements.txt"


    hap_order=read_haplotype_order(order_file)

    hap_edges=parse_unique_elements(unique_file)


    updater=GraphUpdater(base_gfa)

    tool=TOOLS[tool_name]()

    out_dir = repo / "output" / dataset_path / tool.name / MODE_DIRS[mode]

    out_dir.mkdir(parents=True, exist_ok=True)

    time_path = out_dir / f"{tool.name}_time.txt"

    sb_path = out_dir / f"{tool.name}_superbubbles.txt"

    try:

        with open(time_path,"w") as timefile, open(sb_path,"w") as sbfile:

            run_updates(
                mode,
                updater,
                hap_order,
                hap_edges,
                tool,
                timefile,
                sbfile
            )

    finally:

        updater.cleanup()

        tool.close()

    print("Times written to:",time_path)

    print("Superbubbles written to:",sb_path)


def datasets_to_run(repo, d, chromosomes):
    """Dataset paths (relative to input/) for -d and -c; Human without -c = every prepared chromosome."""

    paths = []

    if d in ("Ecoli", "all"):
        paths.append(pathlib.Path("Ecoli"))

    # Human and HPRC (Minigraph-Cactus graphs prepared by mcToDag.py): one
    # folder per chromosome; "all" = Ecoli + Human (HPRC is run on its own)
    for group in (["Human"] if d == "all" else [d] if d in ("Human", "HPRC") else []):
        if chromosomes:
            names = [c if c.startswith("chromosome_") else f"chromosome_{c}" for c in chromosomes]
        else:
            names = sorted(p.name for p in (repo / "input" / group).glob("chromosome_*") if p.is_dir())
        skipped = []
        for name in names:
            path = pathlib.Path(group) / name
            if (repo / "input" / path / "Graphs" / "reduced_graph.gfa").exists():
                paths.append(path)
            else:
                skipped.append(name)
        if skipped:
            print(f"[!] skipped, no reduced_graph.gfa (not preprocessed): {' '.join(skipped)}")

    return paths


#######################################################

if __name__=="__main__":

    parser=argparse.ArgumentParser(
        description="Run a superbubble tool after every update of a dataset; "
                    "see the guide at the top of this file.")

    parser.add_argument(
        "-d",
        choices=["Ecoli","Human","HPRC","all"],
        required=True,
        help="Dataset: Ecoli, Human, HPRC (Minigraph-Cactus, experimental), "
             "or all (= Ecoli + every human chromosome)"
    )

    parser.add_argument("-c",
        nargs="+",
        help="Human only: chromosomes, e.g. 14 X or chromosome_14 "
             "(default: every chromosome that has been preprocessed)")

    parser.add_argument("-m",
        choices=["perEdge","perBatch"],
        required=True,
        help="perEdge: one tool run per edge; perBatch: one per haplotype"
    )

    parser.add_argument(
        "-t",
        "--tool",
        choices=sorted(TOOLS),
        default="bubblegun",
        help="Superbubble tool to run; also names the output directory"
    )

    args=parser.parse_args()

    repo = REPO_ROOT_PATH

    datasets = datasets_to_run(repo, args.d, args.c)

    # One dataset after the other; a dataset that fails (e.g. the tool
    # reports a warning) is reported and the others still run.
    status = {}

    for dataset_path in datasets:

        print(f"\n========== {dataset_path} ({args.m}) ==========", flush=True)

        start = time.time()

        try:
            run_dataset(repo, dataset_path, args.m, args.tool)
            status[dataset_path] = f"ok ({time.time() - start:.0f} s)"
        except (Exception, SystemExit) as e:
            status[dataset_path] = f"FAILED: {e}"
            print(f"[!] {dataset_path} FAILED: {e}", flush=True)

    print("\n===================\nSUMMARY\n===================")

    for dataset_path, st in status.items():
        print(f"{str(dataset_path):22s} {st}")

    sys.exit(1 if any(st.startswith("FAILED") for st in status.values()) else 0)
