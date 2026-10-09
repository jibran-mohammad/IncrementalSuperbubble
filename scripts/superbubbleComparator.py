"""
superbubbleComparator.py

Runs an external superbubble tool (BubbleGun, ...) on the same
graph updates that the C++ algorithms process, for comparison.

Features:

✔ Dynamic repository discovery
✔ Dataset aware
✔ Human chromosome support
✔ Per-edge mode
✔ Per-batch mode
✔ Tool-independent architecture (pick the tool with -t)

Datasets:

Ecoli:

input/Ecoli/Graphs/

Human:

input/Human/chromosome_X/Graphs/

Required files:

reduced_graph.gfa
ordered_haplotypes.txt
unique_elements.txt

Output (directory named after the tool):

output/Ecoli/<tool>/<perEdge|perbatch>/
output/Human/chromosome_X/<tool>/<perEdge|perbatch>/

    <tool>_time.txt          one row per update (edge or haplotype):
                             time of one tool run on the whole graph,
                             its phases (BubbleGun: read, find, chains),
                             and the number of bubbles of each type,
                             followed by the total times
    <tool>_superbubbles.txt  every bubble the tool reports, for the
                             base graph and after each haplotype
                             (the last block is the final graph)


Usage Examples:

Ecoli per-edge:

python superbubbleComparator.py -d Ecoli -m perEdge


Human chromosome 14 per-batch:

python superbubbleComparator.py -d Human -c 14 -m perBatch -t bubblegun
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

def write_time_row(timefile,step,hap,update,result):

    counts = result["counts"]

    phases = result.get("phases", {})

    timefile.write(
        f"{step},{hap},{update},{result['time']},"
        + "".join(f"{t}," for t in phases.values())
        + ",".join(str(n) for n in counts.values())
        + f",{sum(counts.values())}\n"
    )


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

    result = tool.run(updater.write_temp_gfa(), keep_bubbles=True)

    base_time = result["time"]

    timefile.write(
        "Step,Haplotype," + ("Edge" if mode == "perEdge" else "Edges")
        + ",Time,"
        + "".join(f"{p.capitalize()}Time," for p in result.get("phases", {}))
        + ",".join(result["counts"]) + ",TotalBubbles\n"
    )

    write_time_row(timefile,"Base","-","-",result)

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

                result = tool.run(updater.write_temp_gfa(), keep_bubbles=last)

                update_time += result["time"]

                write_time_row(timefile,step,hap,edge,result)

        else:

            for edge in edges:
                updater.add_edge(edge)

            step += 1

            result = tool.run(updater.write_temp_gfa(), keep_bubbles=True)

            update_time += result["time"]

            write_time_row(timefile,step,hap,len(edges),result)

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


#######################################################

if __name__=="__main__":

    parser=argparse.ArgumentParser()

    parser.add_argument(
        "-d",
        choices=["Ecoli","Human"],
        required=True
    )

    parser.add_argument("-c")

    parser.add_argument("-m",
        choices=["perEdge","perBatch"],
        required=True
    )

    parser.add_argument(
        "-t",
        "--tool",
        choices=sorted(TOOLS),
        default="bubblegun",
        help="Superbubble tool to run; also names the output directory"
    )

    args=parser.parse_args()


    ###################################################

    repo = REPO_ROOT_PATH


    ###################################################

    if args.d == "Ecoli":

        dataset_path = pathlib.Path("Ecoli")

    elif args.d == "Human":

        if args.c is None:

            print("❌ Error: Human dataset requires chromosome argument.")
            print("Example:")
            print("python scripts/superbubbleComparator.py -d Human -c chromosome_14 -m perEdge")

            sys.exit(1)

        args.c = f"chromosome_{args.c}" if not args.c.startswith("chromosome_") else args.c    

        dataset_path = pathlib.Path("Human") / args.c

    graph_dir = repo / "input" / dataset_path / "Graphs"


    ###################################################

    base_gfa=graph_dir/"reduced_graph.gfa"

    order_file=graph_dir/"ordered_haplotypes.txt"

    unique_file=graph_dir/"unique_elements.txt"


    hap_order=read_haplotype_order(order_file)

    hap_edges=parse_unique_elements(unique_file)


    updater=GraphUpdater(base_gfa)

    tool=TOOLS[args.tool]()

    out_dir = repo / "output" / dataset_path / tool.name / MODE_DIRS[args.m]

    out_dir.mkdir(parents=True, exist_ok=True)

    time_path = out_dir / f"{tool.name}_time.txt"

    sb_path = out_dir / f"{tool.name}_superbubbles.txt"

    try:

        with open(time_path,"w") as timefile, open(sb_path,"w") as sbfile:

            run_updates(
                args.m,
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
