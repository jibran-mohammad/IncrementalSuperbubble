"""
pipelineRunner.py

Measures the end-to-end time of the best algorithm on a dataset: from the
pangenome graph (GFA) through preprocessing to the finished superbubble
computation. This is the time used to compare the algorithm with
external tools (Preprocessing/toolComparisonToLatex.py).

=====================================================================
HOW TO RUN
=====================================================================

Run from the repository root, with the virtual environment active
(source venv/bin/activate). The graph must exist first
(scripts/graphRunner.py builds it).

    # best algorithm (BEST_ALGORITHM below) on E. coli
    python3 scripts/pipelineRunner.py -d Ecoli -m perEdge

    # per haplotype instead of per edge
    python3 scripts/pipelineRunner.py -d Ecoli -m perBatch

    # human chromosomes 1 and X (note the format: c1, cX)
    python3 scripts/pipelineRunner.py -d Human -c c1 cX -m perEdge

    # all human chromosomes
    python3 scripts/pipelineRunner.py -d Human -m perEdge

    # other algorithms, more timed runs
    python3 scripts/pipelineRunner.py -d Ecoli -m perEdge -a static incBST --repeat 20

Options:

    -d           Ecoli | Human                                   (default: Ecoli)
    -a           algorithm IDs, see ALGORITHMS.md        (default: BEST_ALGORITHM)
    -c           Human only: c1 ... c22, cX, cY, cM      (default: all chromosome_*
                                                         folders in input/Human/;
                                                         those without a graph are
                                                         skipped with a warning)
    -m           perEdge | perBatch                              (default: perBatch)
    --repeat     timed runs per algorithm; AlgorithmTime is their mean (default: 10)

To change the algorithm that runs by default, edit BEST_ALGORITHM below.
Algorithms in NOT_READY_ALGORITHMS (incBSTBatch) are never run.

=====================================================================
WHAT IT DOES
=====================================================================

For each graph (E. coli: the original GFA in input/Ecoli/Graphs/, not
_mg.gfa, _pos.gfa or reduced_graph.gfa; Human:
input/Human/chromosome_<X>/Graphs/chromosome_<X>.gfa):

  1. Preprocessing, timed once: normalizeOrientations.py,
     haplotypeRemoval.py, haplotypeSorter.py, gfaToEdgeListConverter.py,
     edgeListProcessor.py (GFA -> numeric_edgelist.txt). This rewrites
     the files in the dataset's Graphs/ folder (with the same content if
     the graph has not changed).
  2. For each algorithm: superbubbleRunner.py in pipeline mode compiles it
     with BENCHMARK=1 (no per-step timing, no debug output, no output
     files) and runs it --repeat times. Each time is the wall-clock time
     of the program itself: reading the edge list, building the base
     graph and processing all updates. Compilation is not timed.

=====================================================================
OUTPUT
=====================================================================

    output/Ecoli/incsuperbubble/<perEdge|perbatch>/<ID>_time.txt
    output/Human/chromosome_<X>/incsuperbubble/<perEdge|perbatch>/<ID>_time.txt

        PreprocessingTime   GFA -> edge list (step 1)
        Runs                number of timed runs (--repeat)
        AlgorithmTime       mean time of one run (step 2)
        AlgorithmTimeStd    standard deviation of the runs
        Time                PreprocessingTime + AlgorithmTime

"incsuperbubble" (TOOL_NAME below) is the folder of this repository's
algorithms, next to the folders of external tools (e.g. bubblegun/).
Peak memory and the superbubble list are not measured here (BENCHMARK
build); they come from scripts/superbubbleRunner.py.
"""

import os
import sys
import time
import argparse
import subprocess
import pathlib
import statistics
import tempfile
# =============================
# Dynamic repo discovery
# =============================

def find_repo_root(start_path, repo_name="IncSuperbubble"):

    current = start_path.resolve()

    while current != current.parent:
        if current.name == repo_name:
            return current
        current = current.parent

    raise RuntimeError("Repository root not found")


SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
REPO_ROOT = find_repo_root(SCRIPT_DIR)

PREPROCESS_DIR = REPO_ROOT / "Preprocessing"
SCRIPT_DIR = REPO_ROOT / "scripts"

INPUT_DIR = REPO_ROOT / "input"
BASE_OUTPUT_DIR = REPO_ROOT / "output"

# Results go to output/<dataset>[/<chromosome>]/<TOOL_NAME>/<mode>/,
# next to the folders of external tools such as output/Ecoli/bubblegun/.
TOOL_NAME = "incsuperbubble"

# Algorithm run by default (an ID from ALGORITHMS in superbubbleRunner.py;
# see ALGORITHMS.md). incMinIBST = algorithms/incMinIBST.cpp
BEST_ALGORITHM = "incMinIBST"

# Algorithms that are never run by the pipeline (not ready yet).
NOT_READY_ALGORITHMS = {"incBSTBatch"}

# Same mode directory names as superbubbleRunner.py
MODE_DIRS = {"perEdge": "perEdge", "perBatch": "perbatch"}

SUPERBUBBLE_RUNNER = SCRIPT_DIR / "superbubbleRunner.py"


# =============================
# Preprocessing Pipeline
# =============================

def detect_reference_fasta(dataset):

    if dataset == "Ecoli":

        fasta_dir = INPUT_DIR / "Ecoli" / "coli3682_dataset"

    elif dataset == "Human":

        # Reference fasta inside chromosome directory
        return None


    for f in os.listdir(fasta_dir):

        if f.startswith("GCF_") and (
            f.endswith(".fna")
            or f.endswith(".fa")
            or f.endswith(".fasta")
        ):

            return f


    raise RuntimeError(
        "Reference FASTA not found in " + str(fasta_dir)
    )

def run_preprocessing(gfa_file, dataset):

    normalize = PREPROCESS_DIR / "normalizeOrientations.py"
    remove = PREPROCESS_DIR / "haplotypeRemoval.py"
    sorter = PREPROCESS_DIR / "haplotypeSorter.py"
    converter = PREPROCESS_DIR / "gfaToEdgeListConverter.py"
    processor = PREPROCESS_DIR / "edgeListProcessor.py"


    # ---------- Environment variables ----------
    env = os.environ.copy()

    if dataset == "Human":
        chrom = pathlib.Path(gfa_file).parent.parent.name
        env["PIPELINE_DATASET_NAME"] = f"Human/{chrom}"
    else:
        env["PIPELINE_DATASET_NAME"] = dataset


    # Derive reference fasta from GFA name
    gfa_name = os.path.basename(gfa_file)

    base = gfa_name.replace(".gfa","")

    ref_fasta = base + ".fna"

    env["PIPELINE_REF_FASTA"] = ref_fasta


    print("\nUsing reference fasta:")
    print(ref_fasta)
    print()


    cmds = [

        ["python3", normalize, gfa_file],
        ["python3", remove],
        ["python3", sorter],
        ["python3", converter],
        ["python3", processor]

    ]

    for cmd in cmds:

        print("Running:", " ".join(map(str,cmd)))

        subprocess.run(cmd, check=True, env=env)



# =============================
# Run superbubble
# =============================

def run_superbubble(dataset,algorithm,mode,chromosomes,repeat):
    """
    Runs one algorithm through superbubbleRunner.py in pipeline mode
    (BENCHMARK build) `repeat` times and returns the wall-clock times of
    the executable itself (one per run), i.e. without make clean /
    compilation. Each time covers reading the edge list, building the base
    graph and processing all updates.
    """

    cmd = [

        "python3",
        SUPERBUBBLE_RUNNER,
        "-d",dataset,
        "-m",mode,
        "-a",algorithm

    ]

    if chromosomes:
        cmd += ["-c"]+chromosomes


    fd, timing_file = tempfile.mkstemp(prefix="pipeline_timing_", suffix=".txt")
    os.close(fd)

    env = os.environ.copy()
    env["SUPERBUBBLE_PIPELINE_MODE"] = "1"
    env["SUPERBUBBLE_TIMING_FILE"] = timing_file
    env["SUPERBUBBLE_PIPELINE_RUNS"] = str(repeat)

    try:

        subprocess.run(cmd, check=True, env=env)

        with open(timing_file) as f:
            times = [float(line.split("\t")[1]) for line in f if line.strip()]

    finally:

        os.remove(timing_file)

    if len(times) != repeat:
        raise RuntimeError(f"Expected {repeat} execution times for {algorithm}, got {len(times)}")

    return times



# =============================
# Dataset discovery
# =============================

def discover_gfa_files(dataset, chromosomes):
    """
    Returns a list of (gfa_file, chromosome_name or None).
    Only original graphs are used, never intermediate GFAs.
    """

    gfas = []

    forbidden_suffixes = [
        "_mg.gfa",
        "_pos.gfa",
        "reduced_graph.gfa",
        "_normalized.gfa",
        "_sorted.gfa"
    ]

    if dataset == "Ecoli":

        gfa_dir = INPUT_DIR / "Ecoli" / "Graphs"

        for f in sorted(os.listdir(gfa_dir)):

            if f.endswith(".gfa") and not any(f.endswith(x) for x in forbidden_suffixes):

                gfas.append((gfa_dir / f, None))

    elif dataset == "Human":

        base_dir = INPUT_DIR / "Human"

        if chromosomes:
            chrs = chromosomes
        else:
            # AUTO: pick all chromosomes
            chrs = sorted(
                d for d in os.listdir(base_dir)
                if d.startswith("chromosome_")
            )

        for c in chrs:

            gfa_dir = base_dir / c / "Graphs"
            expected = f"{c}.gfa"
            gfa_file = gfa_dir / expected

            # Chromosomes without a graph are skipped (like graphRunner.py and
            # superbubbleRunner.py do), so one missing folder does not stop
            # the whole run.
            if not gfa_file.exists():
                print(f"[!] Warning: no graph for {c} ({gfa_file}), skipping it.", file=sys.stderr)
                continue

            gfas.append((gfa_file, c))

    if len(gfas) == 0:

        raise RuntimeError(
            "No valid input .gfa file found.\n"
            "Expected original .gfa (not _mg.gfa, _pos.gfa, reduced_graph.gfa)"
        )

    if dataset == "Ecoli" and len(gfas) > 1:

        print("\nWARNING: Multiple candidate GFAs found, using the first one:")

        for g, _ in gfas:
            print("  ", g)

        gfas = gfas[:1]

    print("\nSelected input GFA:")
    for g, _ in gfas:
        print("  ", g)
    print()

    return gfas

def normalize_chromosome_names(chromosomes):

    if not chromosomes:
        return None

    normalized=[]

    for c in chromosomes:

        c=c.strip()

        if not c.lower().startswith("c"):
            raise ValueError(
                f"Invalid chromosome format: {c}\n"
                "Use c1,c2,...,c22,cX,cY,cM"
            )

        name=c[1:]

        if name.isdigit():

            normalized.append(f"chromosome_{name}")

        elif name.upper() in ["X","Y","M"]:

            normalized.append(f"chromosome_{name.upper()}")

        else:

            raise ValueError(
                f"Invalid chromosome: {c}"
            )

    return normalized

def build_output_dir(dataset, chrom, mode):

    out_dir = BASE_OUTPUT_DIR / dataset

    if chrom:
        out_dir = out_dir / chrom

    return out_dir / TOOL_NAME / MODE_DIRS[mode]

# =============================
# Main
# =============================

def main():

    parser=argparse.ArgumentParser(
        description="End-to-end time of the best algorithm: GFA -> preprocessing -> "
                    "superbubbles. See the guide at the top of this file.")

    parser.add_argument("-d",default="Ecoli",choices=["Ecoli","Human"],
                        help="Dataset (default: Ecoli)")

    parser.add_argument("-a", nargs="+",
                        help=f"Algorithms to run (default: {BEST_ALGORITHM})")

    parser.add_argument("-c",nargs="+",
                        help="Human chromosomes: c1 ... c22, cX, cY, cM (default: all)")

    parser.add_argument("-m",default="perBatch",
                        choices=["perBatch","perEdge"],
                        help="perEdge: update after every edge; perBatch: after every haplotype (default: perBatch)")

    parser.add_argument("--repeat", type=int, default=10,
                        help="Timed runs per algorithm; AlgorithmTime is their mean (default: 10)")

    args=parser.parse_args()


    from superbubbleRunner import ALGORITHMS

    AVAILABLE_ALGORITHMS = [algo["name"] for algo in ALGORITHMS]

    algorithms = args.a if args.a else [BEST_ALGORITHM]

    for algo in algorithms:
        if algo not in AVAILABLE_ALGORITHMS:
            raise ValueError(f"Invalid algorithm: {algo}")

    skipped = [algo for algo in algorithms if algo in NOT_READY_ALGORITHMS]
    algorithms = [algo for algo in algorithms if algo not in NOT_READY_ALGORITHMS]

    for algo in skipped:
        print(f"[SKIP] {algo} is not ready yet")

    if not algorithms:
        raise ValueError("No algorithm left to run")

    chromosomes=normalize_chromosome_names(args.c)

    gfa_files=discover_gfa_files(args.d,chromosomes)

    total_start=time.time()


    for gfa, chrom_name in gfa_files:

        # ---------- GFA -> edge list (shared by all algorithms) ----------
        print("\nProcessing:",gfa)

        start = time.time()

        run_preprocessing(str(gfa), args.d)

        preprocessing_time = time.time() - start


        # ---------- Each algorithm on the edge list ----------
        for algo in algorithms:

            times = run_superbubble(
                args.d,
                algo,
                args.m,
                [chrom_name] if chrom_name else None,
                args.repeat
            )

            execution_time = statistics.mean(times)

            execution_std = statistics.stdev(times) if len(times) > 1 else 0.0

            out_dir = build_output_dir(args.d, chrom_name, args.m)

            out_dir.mkdir(parents=True, exist_ok=True)

            out_file = out_dir / f"{algo}_time.txt"

            with open(out_file, "w") as f:
                f.write(f"Dataset: {args.d}\n")
                if chrom_name:
                    f.write(f"Chromosome: {chrom_name}\n")
                f.write(f"Algorithm: {algo}\n")
                f.write(f"Mode: {args.m}\n")
                f.write(f"GFA: {gfa}\n")
                f.write(f"PreprocessingTime: {preprocessing_time:.6f} seconds\n")
                f.write(f"Runs: {len(times)}\n")
                f.write(f"AlgorithmTime: {execution_time:.6f} seconds\n")
                f.write(f"AlgorithmTimeStd: {execution_std:.6f} seconds\n")
                f.write(f"Time: {preprocessing_time + execution_time:.6f} seconds\n")

            print(f"\n{algo}: preprocessing {preprocessing_time:.3f}s + algorithm {execution_time:.3f}s"
                  f" = {preprocessing_time + execution_time:.3f}s  -> {out_file}")


    total_time=time.time()-total_start


    print("\n===================")
    print("TOTAL PIPELINE TIME")
    print("===================")

    print(total_time,"seconds")


if __name__=="__main__":
    main()
