"""
superbubbleComparator.py

Runner for BubbleGun comparison experiments.

Features:

✔ Dynamic repository discovery
✔ Dataset aware
✔ Human chromosome support
✔ Per-edge mode
✔ Per-batch mode
✔ Tool-independent architecture

Datasets:

Ecoli:

input/Ecoli/Graphs/

Human:

input/Human/chromosome_X/Graphs/

Required files:

reduced_graph.gfa
ordered_haplotypes.txt
unique_elements.txt


Usage Examples:

Ecoli per-edge:

python superbubbleComparator.py -d Ecoli -m perEdge -o ecoli_perEdge.csv


Human chromosome 14 per-batch:

python superbubbleComparator.py -d Human -c 14 -m perBatch -o ecoli_perEdge.csv
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


#######################################################

def run_perEdge(updater,hap_order,hap_edges,tool,outfile):

    total_runtime = 0

    ################################################
    # Base graph
    ################################################

    gfa = updater.write_temp_gfa()

    runtime, sb = tool.run(gfa)

    total_runtime += runtime

    outfile.write(f"Base,-,-,{runtime},{sb}\n")

    print("Base completed")


    ################################################
    # Edge updates
    ################################################

    for hap in hap_order:

        for edge in hap_edges[hap]:

            updater.add_edge(edge)

            gfa = updater.write_temp_gfa()

            runtime, sb = tool.run(gfa)

            total_runtime += runtime

            outfile.write(
                f"Edge,{hap},{edge},{runtime},{sb}\n"
            )


    ################################################
    # Total BubbleGun runtime
    ################################################

    outfile.write("\n")

    outfile.write(
        f"TotalBubbleGunTime,{total_runtime}\n"
    )

    print("Total BubbleGun time:", total_runtime)


#######################################################

def run_perBatch(updater,hap_order,hap_edges,tool,outfile):

    total_runtime = 0

    ################################################
    # Base graph
    ################################################

    gfa = updater.write_temp_gfa()

    runtime, sb = tool.run(gfa)

    total_runtime += runtime

    outfile.write(f"Base,-,-,{runtime},{sb}\n")

    print("Base completed")


    ################################################
    # Batch updates
    ################################################

    for hap in hap_order:

        for edge in hap_edges[hap]:

            updater.add_edge(edge)

        gfa = updater.write_temp_gfa()

        runtime, sb = tool.run(gfa)

        total_runtime += runtime

        outfile.write(
            f"Batch,{hap},-,{runtime},{sb}\n"
        )


    ################################################
    # Total BubbleGun runtime
    ################################################

    outfile.write("\n")

    outfile.write(
        f"TotalBubbleGunTime,{total_runtime}\n"
    )

    print("Total BubbleGun time:", total_runtime)


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
        "-o",
        "--output",
        default="superbubble_results.csv",
        help="Output file for results"
    )

    args=parser.parse_args()


    ###################################################

    script_dir=pathlib.Path(__file__).parent

    repo = REPO_ROOT_PATH


    ###################################################

    if args.d == "Ecoli":

        graph_dir = repo / "input" / "Ecoli" / "Graphs"

    elif args.d == "Human":

        if args.c is None:

            print("❌ Error: Human dataset requires chromosome argument.")
            print("Example:")
            print("python scripts/superbubbleComparator.py -d Human -c chromosome_14 -m perEdge")

            sys.exit(1)

        args.c = f"chromosome_{args.c}" if not args.c.startswith("chromosome_") else args.c    

        graph_dir = repo / "input" / "Human" / args.c / "Graphs"


    ###################################################

    base_gfa=graph_dir/"reduced_graph.gfa"

    order_file=graph_dir/"ordered_haplotypes.txt"

    unique_file=graph_dir/"unique_elements.txt"


    hap_order=read_haplotype_order(order_file)

    hap_edges=parse_unique_elements(unique_file)


    updater=GraphUpdater(base_gfa)

    tool=BubbleGunTool()

    output_file = args.output

    with open(output_file,"w") as outfile:

        outfile.write(
            "Mode,Haplotype,Edge,Runtime,Superbubbles\n"
        )

        if args.m=="perEdge":

            run_perEdge(
                updater,
                hap_order,
                hap_edges,
                tool,
                outfile
            )

        else:

            run_perBatch(
                updater,
                hap_order,
                hap_edges,
                tool,
                outfile
            )

    print("Results written to:",output_file)