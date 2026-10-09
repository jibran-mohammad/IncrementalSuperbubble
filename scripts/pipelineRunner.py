import os
import sys
import time
import argparse
import subprocess
import pathlib
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
OUTPUT_DIR = BASE_OUTPUT_DIR / "pipelineRunner"

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

def run_superbubble(dataset,algorithm,mode,chromosomes):
    """
    Runs one algorithm through superbubbleRunner.py in pipeline mode
    (BENCHMARK build, single run) and returns the wall-clock time of
    the executable itself, i.e. without make clean / compilation.
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

    try:

        subprocess.run(cmd, check=True, env=env)

        with open(timing_file) as f:
            times = [float(line.split("\t")[1]) for line in f if line.strip()]

    finally:

        os.remove(timing_file)

    if not times:
        raise RuntimeError(f"No execution time reported for {algorithm}")

    return sum(times)



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

            if not gfa_file.exists():
                raise RuntimeError(f"Missing GFA: {gfa_file}")

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

def build_output_filename(dataset, chrom, algorithm, mode):

    parts = [dataset]

    if chrom:
        parts.append(chrom)

    if algorithm:
        parts.append(algorithm)

    parts.append(mode)

    return "_".join(parts) + ".txt"

# =============================
# Main
# =============================

def main():

    parser=argparse.ArgumentParser()

    parser.add_argument("-d",default="Ecoli",choices=["Ecoli","Human"])

    parser.add_argument("-a", nargs="+")

    parser.add_argument("-c",nargs="+")

    parser.add_argument("-m",default="perBatch",
                        choices=["perBatch","perEdge"])

    args=parser.parse_args()


    from superbubbleRunner import ALGORITHMS

    AVAILABLE_ALGORITHMS = [algo["name"] for algo in ALGORITHMS]

    if args.a:
        algorithms = args.a
        for algo in algorithms:
            if algo not in AVAILABLE_ALGORITHMS:
                raise ValueError(f"Invalid algorithm: {algo}")
    else:
        algorithms = AVAILABLE_ALGORITHMS


    os.makedirs(OUTPUT_DIR, exist_ok=True)
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

            if args.m == "perEdge" and algo == "IncrementalSuperbubbleBatch":
                print(f"[SKIP] {algo} not valid for perEdge")
                continue

            execution_time = run_superbubble(
                args.d,
                algo,
                args.m,
                [chrom_name] if chrom_name else None
            )

            filename = build_output_filename(
                args.d,
                chrom_name,
                algo,
                args.m
            )

            out_file = OUTPUT_DIR / filename

            with open(out_file, "w") as f:
                f.write(f"Dataset: {args.d}\n")
                if chrom_name:
                    f.write(f"Chromosome: {chrom_name}\n")
                f.write(f"Algorithm: {algo}\n")
                f.write(f"Mode: {args.m}\n")
                f.write(f"GFA: {gfa}\n")
                f.write(f"PreprocessingTime: {preprocessing_time:.6f} seconds\n")
                f.write(f"AlgorithmTime: {execution_time:.6f} seconds\n")
                f.write(f"Time: {preprocessing_time + execution_time:.6f} seconds\n")

            print(f"\n{algo}: preprocessing {preprocessing_time:.3f}s + algorithm {execution_time:.3f}s"
                  f" = {preprocessing_time + execution_time:.3f}s  -> {out_file}")


    total_time=time.time()-total_start


    print("\n===================")
    print("TOTAL PIPELINE TIME")
    print("===================")

    print(total_time,"seconds")


    total_file = OUTPUT_DIR / "TOTAL_PIPELINE_TIME.txt"

    with open(total_file, "w") as f:
        f.write(f"Total Pipeline Time: {total_time:.6f} seconds\n")


if __name__=="__main__":
    main()


"""
Times written per (dataset, chromosome, algorithm, mode) to output/pipelineRunner/:
  PreprocessingTime  GFA -> edge list (normalize, haplotype removal/sort, conversion)
  AlgorithmTime      one run of the compiled algorithm (BENCHMARK build, compilation excluded)
  Time               PreprocessingTime + AlgorithmTime

python3 pipelineRunner.py -d Human -a SimpleSuperbubble -m perEdge
Runs all chromosomes for SimpleSuperbubble
python3 pipelineRunner.py -d Ecoli -a SimpleSuperbubble -m perEdge
python3 pipelineRunner.py -d Human -c c1 -a SimpleSuperbubble -m perEdge


"""
