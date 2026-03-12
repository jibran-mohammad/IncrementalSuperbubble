import os
import sys
import time
import argparse
import subprocess
import pathlib

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
OUTPUT_DIR = REPO_ROOT / "output"

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

    cmd = [

        "python3",
        SUPERBUBBLE_RUNNER,
        "-d",dataset,
        "-m",mode

    ]

    if algorithm:
        cmd += ["-a",algorithm]

    if chromosomes:
        cmd += ["-c"]+chromosomes


    env = os.environ.copy()
    env["SUPERBUBBLE_PIPELINE_MODE"] = "1"

    subprocess.run(cmd, check=True, env=env)



# =============================
# Dataset discovery
# =============================

def discover_gfa_files(dataset, chromosomes):

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

        for f in os.listdir(gfa_dir):

            if f.endswith(".gfa") and not any(f.endswith(x) for x in forbidden_suffixes):

                gfas.append(gfa_dir / f)

    elif dataset == "Human":

        if not chromosomes:
            raise RuntimeError(
                "Human dataset requires chromosome selection (-c c1,c2,...,cX)"
            )

        for c in chromosomes:

            gfa_dir = INPUT_DIR / "Human" / c / "Graphs"
            expected = f"{c}.gfa"

            gfa_file = gfa_dir / expected

            if not gfa_file.exists():
                raise RuntimeError(
                    f"Expected GFA file not found:\n{gfa_file}"
                )

            # Safety check (paranoia mode)
            if any(expected.endswith(x) for x in forbidden_suffixes):
                raise RuntimeError(
                    f"Invalid GFA selected (intermediate file): {gfa_file}"
                )

            gfas.append(gfa_file)

    print("\nSelected input GFA:")
    for g in gfas:
        print("  ", g)
    print()

    return gfas

    gfas = []


    if dataset == "Ecoli":

        gfa_dir = INPUT_DIR / "Ecoli" / "Graphs"

        for f in os.listdir(gfa_dir):

            if is_original_gfa(f):

                gfas.append(gfa_dir / f)


    elif dataset == "Human":

        if chromosomes:
            chrs = chromosomes
        else:
            chrs = [d for d in os.listdir(INPUT_DIR / "Human") if d.startswith("chromosome_")]


        for c in chrs:

            gfa_dir = INPUT_DIR / "Human" / c / "Graphs"

            for f in os.listdir(gfa_dir):

                if is_original_gfa(f):

                    gfas.append(gfa_dir / f)


    if len(gfas) == 0:

        raise RuntimeError(
            "No valid input .gfa file found.\n"
            "Expected original .gfa (not _mg.gfa, _pos.gfa, reduced_graph.gfa)"
        )


    if len(gfas) > 1:

        print("\nWARNING: Multiple candidate GFAs found:")

        for g in gfas:
            print("  ", g)

        print("\nUsing first one.\n")

    print("\nSelected input GFA:")
    for g in gfas:
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

# =============================
# Main
# =============================

def main():

    parser=argparse.ArgumentParser()

    parser.add_argument("-d",default="Ecoli")

    parser.add_argument("-a")

    parser.add_argument("-c",nargs="+")

    parser.add_argument("-m",default="perBatch",
                        choices=["perBatch","perEdge"])

    args=parser.parse_args()


    chromosomes=normalize_chromosome_names(args.c)

    gfa_files=discover_gfa_files(args.d,chromosomes)

    total_start=time.time()


    for gfa in gfa_files:

        print("\nProcessing:",gfa)

        run_preprocessing(str(gfa), args.d)


    run_superbubble(
        args.d,
        args.a,
        args.m,
        chromosomes
    )

    total_time=time.time()-total_start


    print("\n===================")
    print("TOTAL PIPELINE TIME")
    print("===================")

    print(total_time,"seconds")


    out_file=OUTPUT_DIR/"TotalPipelineTime.txt"

    with open(out_file,"w") as f:

        f.write(str(total_time))


if __name__=="__main__":
    main()


"""
python3 pipelineRunner.py -d Ecoli -a SimpleSuperbubble -m perBatch
python3 pipelineRunner.py -d Ecoli -a SimpleSuperbubble -m perEdge
python3 pipelineRunner.py -d Human -c c1 -a SimpleSuperbubble -m perEdge


"""    