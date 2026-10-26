"""
superbubbleRunner.py

Compiles the C++ superbubble algorithms and runs them on the edge list of
a dataset, 100 times each; writes the averaged timing file and the final
superbubble list of every algorithm. These files feed the plots and
tables (parse_to_csv.py, plotcsv.py, csvToLatex.py) and the comparison
with external tools (toolComparisonToLatex.py).

=====================================================================
HOW TO RUN
=====================================================================

Run from the repository root, with the virtual environment active
(source venv/bin/activate). The edge list must exist first
(scripts/graphRunner.py).

    # E. coli, all algorithms, one update per edge
    python3 scripts/superbubbleRunner.py -d Ecoli -m perEdge

    # E. coli, one update per haplotype
    python3 scripts/superbubbleRunner.py -d Ecoli -m perBatch

    # only some algorithms (IDs: see ALGORITHMS.md)
    python3 scripts/superbubbleRunner.py -d Ecoli -m perEdge -a static incMinIBST

    # all human chromosomes: chromosome_1 ... chromosome_22, chromosome_X,
    # chromosome_Y, chromosome_M; those without an edge list are skipped
    # with a warning
    python3 scripts/superbubbleRunner.py -d Human -m perEdge

    # specific human chromosomes (full folder names)
    python3 scripts/superbubbleRunner.py -d Human -m perEdge -c chromosome_1 chromosome_X

Options:

    --runs N             runs per algorithm, averaged             (default: 100;
                         1 for very large graphs such as HPRC chromosomes)
    -d / --dataset       Ecoli | Human | HPRC                     (default: Ecoli;
                         HPRC: Minigraph-Cactus graphs prepared by
                         Preprocessing/mcToDag.py, experimental)
    -a / --algorithms    algorithm IDs: static incRecompute incNormalBST incBST incMinIBST
                         incBSTBatch                              (default: all)
    -c / --chromosomes   Human only: chromosome_1 chromosome_X ... (default: all)
    -m / --mode          perEdge: recompute after every edge;
                         perBatch: after every haplotype           (default: perBatch)

Notes:

  - incBSTBatch (the batch algorithm) only runs in perBatch mode (skipped
    in perEdge): it processes a whole haplotype at once.
  - incNormalBST only runs in perEdge mode (skipped in perBatch): it
    re-keys its trees after every edge, so in perBatch it is far slower
    than every other algorithm (Ecoli: 55 s vs 0.4 s static) and adds
    nothing there.
  - Every run starts with "make clean" and uses one _tmp_runs/ folder in
    the repository root, so run only one superbubbleRunner.py at a time
    per copy of the repository.
  - NUM_RUNS (below) sets the number of runs per algorithm (100); the
    times in the output are their average. Debug output is off (build
    with make DEBUG=1 to see it).
  - Run by pipelineRunner.py, this script works in "pipeline mode"
    (SUPERBUBBLE_PIPELINE_MODE=1): BENCHMARK build, no timing files, only
    the wall-clock time of each run is reported back.

=====================================================================
INPUT / OUTPUT
=====================================================================

Input:   input/Ecoli/Graphs/numeric_edgelist.txt
         input/Human/chromosome_<X>/Graphs/numeric_edgelist.txt

Output:  output/Ecoli/timeFile/<perEdge|perbatch>/
         output/Human/chromosome_<X>/timeFile/<perEdge|perbatch>/

    <ID>_time.txt           per update: insert, DFS, topological order,
                            minI, maxO and superbubble times, superbubble
                            count, and (incremental algorithms) cases and
                            affected/explored vertices and edges; last
                            line PeakMemory(KB); averaged over NUM_RUNS
    <ID>_superbubbles.txt   superbubbles of the final graph, one
                            "<entrance> <exit>" per line in input ids

Next step: python3 Preprocessing/parse_to_csv.py --dataset ecoli -m peredge
"""

import os
import resource
import signal
import subprocess
import sys
import shutil
import argparse
import time



# --- DYNAMIC PATH CONFIGURATION ---
# Get the absolute path to this script's directory (e.g., .../IncSuperbubble/scripts)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Get the parent directory (the true project root, e.g., .../IncSuperbubble)
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

# 1. The directory where your Makefile is located.
#    The 'make' command will be run here.
MAKE_DIR = PROJECT_ROOT

PIPELINE_MODE = os.environ.get("SUPERBUBBLE_PIPELINE_MODE") == "1"

# In pipeline mode, pipelineRunner.py passes a file here; the wall-clock
# time of each executable run (compilation excluded) is appended to it.
PIPELINE_TIMING_FILE = os.environ.get("SUPERBUBBLE_TIMING_FILE")

# Number of repetitions for averaging. In pipeline mode pipelineRunner.py
# sets SUPERBUBBLE_PIPELINE_RUNS (its --repeat) and averages the reported times.
NUM_RUNS = int(os.environ.get("SUPERBUBBLE_PIPELINE_RUNS", "1")) if PIPELINE_MODE else 100

print("PIPELINE MODE:", PIPELINE_MODE)

# --- CONFIGURATION ---

# 1. A list describing each of your C++ algorithms.
#    - 'name': A friendly name for logging.
#    - 'target': The TARGET name for the Makefile (e.g., 'simple_superbubble_app').
#    - 'output_filename': The final name for the output file.
#
#    NOTE: The 'output_dir' has been removed from this list because the
#    main() function now generates it dynamically based on the dataset.
#    Every algorithm has one ID, used for the source file
#    (algorithms/<ID>.cpp), the make target and binary (test/<ID>), the
#    name below, the timing file (<ID>_time.txt) and the CSV/table column
#    prefix. See ALGORITHMS.md for what each algorithm does.
ALGORITHMS = [
    {
        "name": "static",
        "target": "static",
        "output_filename": "static_time.txt"
    },

    {
        "name": "incRecompute",
        "target": "incRecompute",
        "output_filename": "incRecompute_time.txt"
    },

    {
        # incBST with normal binary search trees (keyed by position) for
        # minI / maxO: every vertex that moves is re-keyed in its neighbours'
        # trees (lib/DFS/incrementalDFSNormalBST.cpp)
        "name": "incNormalBST",
        "target": "incNormalBST",
        "output_filename": "incNormalBST_time.txt"
    },

   
    #{
    #    "name": "IncrementalAlgorithm",
    #   "target": "incsuperbubble",
    #   "output_filename": "SuperbubbleIncrementalDFSOutput.txt"
    #},
    {
        "name": "incBST",
        "target": "incBST",
        "output_filename": "incBST_time.txt"
    },

    {
        "name": "incMinIBST",
        "target": "incMinIBST",
        "output_filename": "incMinIBST_time.txt"
    },

    {
        # incMinIBST + superbubble heuristic: superbubbles updated locally
        # where possible (lib/findSuperbubble/superbubbleh1.cpp)
        "name": "incMinIBSTH1",
        "target": "incMinIBSTH1",
        "output_filename": "incMinIBSTH1_time.txt"
    },
     
    #{
    #    "name": "SimpleSuperbubbleH1",
    #    "target": "simpleSuperbubbleh1",
    #    "output_filename": "SuperbubbleSimpleTimeOutputH1.txt"
    #},

    {
        "name": "incBSTBatch",
        "target": "incBSTBatch",
        "output_filename": "incBSTBatch_time.txt"
    }
]
# ---------------------

def run_make_clean(cpp_dir):
    """Runs 'make clean' in the specified directory."""
    print(f"-> Running 'make clean' in '{cpp_dir}'...")
    clean_command = ["make", "clean"]
    try:
        # The 'cwd' argument ensures the command runs in the correct directory.
        subprocess.run(
            clean_command,
            cwd=cpp_dir,
            check=True,
            capture_output=True,
            text=True
        )
        print("   ...Previous builds cleaned successfully.")
        return True
    except FileNotFoundError:
        print(f"   [!] ERROR: 'make' command not found. Is it installed and in your system's PATH?", file=sys.stderr)
        return False
    except subprocess.CalledProcessError as e:
        # If 'make clean' fails, it might be because there's nothing to clean,
        # but it's safer to report the error and stop.
        print(f"   [!] 'make clean' FAILED.", file=sys.stderr)
        print("--- make error output ---", file=sys.stderr)
        print(e.stderr, file=sys.stderr)
        print("-------------------------", file=sys.stderr)
        return False

def compile_target_with_make(cpp_dir, target):
    """
    Compiles a target using the Makefile and captures the last line of stdout
    as the path to the executable.
    Returns the executable path on success, None on failure.
    """
    print(f"-> Compiling target '{target}' using Makefile in '{cpp_dir}'...")
    if PIPELINE_MODE:
        make_command = ["make", f"TARGET={target}", "BENCHMARK=1"]
    else:
        make_command = ["make", f"TARGET={target}"]
    try:
        # Run the make command, capturing its output.
        result = subprocess.run(
            make_command,
            cwd=cpp_dir,
            check=True,
            capture_output=True,
            text=True
        )
        # Assume the last line of stdout is the path to the executable.
        # .strip() removes leading/trailing whitespace, including the final newline.
        last_line = result.stdout.strip().splitlines()[-1]
        
        executable_path = None
        if last_line.startswith("Linking "):
            # Split the line by spaces, e.g., "Linking test/path..." -> ["Linking", "test/path..."]
            path_part = last_line.split()[1]
            # Remove the trailing "..." from the path part.
            relative_path = path_part.rstrip('.')
            # Construct the full, absolute path by joining with the base directory.
            executable_path = os.path.join(cpp_dir, relative_path)
        else:
            # Fallback for makefiles that just output the path directly.
            executable_path = last_line
        print(executable_path)
        if not executable_path:
             print(f"   [!] ERROR: Could not parse executable path from make output.", file=sys.stderr)
             return None

        print(f"   ...Success! Parsed executable path: '{executable_path}'")
        return executable_path
    except FileNotFoundError:
        print(f"   [!] ERROR: 'make' command not found. Is it installed and in your system's PATH?", file=sys.stderr)
        return None
    except subprocess.CalledProcessError as e:
        print(f"   [!] Makefile FAILED for target '{target}'.", file=sys.stderr)
        print("--- make error output ---", file=sys.stderr)
        print(e.stderr, file=sys.stderr)
        print("-------------------------", file=sys.stderr)
        return None

def raise_stack_limit():
    """
    Runs in the child process before the C++ program starts: raises its
    stack limit to the largest the system allows (usually unlimited). The
    DFS of the algorithms is recursive, one level per vertex on a path, and
    the default 8 MB stack overflows on very large graphs (segmentation
    fault, e.g. HPRC chromosome 22 with a 740,000-vertex reference path).
    """
    soft, hard = resource.getrlimit(resource.RLIMIT_STACK)
    if soft != hard:
        resource.setrlimit(resource.RLIMIT_STACK, (hard, hard))


def run_cpp_program(executable_path, input_path,full_output_path, mode):
    """Runs a compiled C++ program, redirecting the input file to its stdin."""
    print(f"-> Running {os.path.basename(executable_path)} with input {os.path.basename(input_path)}...in {mode} mode...")
    print(f"   ...Writing output to: {full_output_path}")

    # The executable path is now absolute, so it's the only item in the command.
    run_command = [executable_path, full_output_path, mode]
    try:
        with open(input_path, 'r') as input_file_handle:
            subprocess.run(
                run_command,
                stdin=input_file_handle,
                check=True,
                cwd=PROJECT_ROOT,
                preexec_fn=raise_stack_limit
            )
        print(f"   ...Finished successfully.")
        return True
    except subprocess.CalledProcessError as e:
        # Why it ended: a negative return code means the program was killed
        # by a signal (it prints nothing itself then, e.g. on a segmentation
        # fault); otherwise its exit code. Its own error messages (stderr)
        # are already in the log above, they are not captured here.
        if e.returncode < 0:
            sig = -e.returncode
            try:
                name = signal.Signals(sig).name
            except ValueError:
                name = "unknown signal"
            hints = {"SIGSEGV": "segmentation fault, e.g. a stack overflow or an invalid memory access",
                     "SIGKILL": "killed, e.g. by the system when it ran out of memory (see dmesg)",
                     "SIGABRT": "aborted, e.g. a failed assertion or std::bad_alloc (out of memory)",
                     "SIGBUS": "bus error (invalid memory access)"}
            reason = f"killed by signal {sig} ({name}" + (f": {hints[name]}" if name in hints else "") + ")"
        else:
            reason = f"exit code {e.returncode}"
        print(f"   [!] Execution FAILED for {os.path.basename(executable_path)}: {reason}", flush=True)
        print(f"   [!] Execution FAILED for {os.path.basename(executable_path)}: {reason}", file=sys.stderr, flush=True)
        return False


def main():
    """Main function to orchestrate the compilation and execution pipeline."""
    all_algo_names = [algo['name'] for algo in ALGORITHMS]
    
    # --- Instructions for the user ---
    # (The print statements below also serve as a mini-guide)
    print("--- Superbubble Pipeline Orchestrator ---")
    print(f"Available algorithms: {', '.join(all_algo_names)}")
    
    # --- Set up command-line argument parsing ---
    parser = argparse.ArgumentParser(
        description="Compile and run specified C++ superbubble algorithms.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '-a', '--algorithms',
        nargs='+',  
        choices=all_algo_names,
        metavar='ALGO_NAME',
        help=f"Specify which algorithms to run. Choices: {', '.join(all_algo_names)}. "
             f"\nIf not provided, all algorithms will be run."
    )
    parser.add_argument(
        '-d', '--dataset',
        type=str,
        default='Ecoli',
        choices=['Ecoli', 'Human', 'HPRC'], # Added choices for clarity
        help="Specify which dataset to run (e.g., 'Ecoli', 'Human', 'HPRC'). Default: Ecoli"
    )
    
    # Add a new argument for specifying chromosomes
    parser.add_argument(
        '-c', '--chromosomes',
        nargs='+',
        metavar='CHR_NAME',
        help="Specify which human chromosomes to run (e.g., 'chromosome_1', 'chromosome_X'). \n"
             "Only used if '-d Human' is selected. \n"
             "If not provided, all human chromosomes will be run."
    )
    
    parser.add_argument(
        '-m', '--mode',
        choices=['perBatch', 'perEdge'],
        default='perBatch',
        help="Specify processing mode. 'perBatch' (after haplotype) or 'perEdge' (after every edge)."
    )
    
    parser.add_argument(
        '--runs',
        type=int,
        default=None,
        help="Runs per algorithm, averaged (default: 100). Use 1 for very large graphs "
             "(e.g. HPRC chromosomes), where one run takes hours."
    )
    args = parser.parse_args()

    # --runs overrides the number of runs (not in pipeline mode, where
    # pipelineRunner.py sets it with --repeat)
    global NUM_RUNS
    if args.runs is not None and not PIPELINE_MODE:
        NUM_RUNS = args.runs

    mode_dir_map = {
        "perBatch": "perbatch",
        "perEdge": "perEdge"
    }
    mode_dir = mode_dir_map[args.mode]

    if PIPELINE_MODE:
        print("Running in PIPELINE BENCHMARK mode")
    else:
        print("Running in NORMAL EXPERIMENT mode")

    if args.algorithms:
        algorithms_to_run = [algo for algo in ALGORITHMS if algo['name'] in args.algorithms]
        print(f"--- Running selected algorithms: {[algo['name'] for algo in algorithms_to_run]} ---")
    else:
        algorithms_to_run = ALGORITHMS
        print("--- Starting Superbubble Algorithm Pipeline (all algorithms) ---")

    if args.mode == "perEdge":
        removed = [algo['name'] for algo in algorithms_to_run if algo['name'] == "incBSTBatch"]
        algorithms_to_run = [
            algo for algo in algorithms_to_run
            if algo['name'] != "incBSTBatch"
        ]
        if removed:
            print("[INFO] Skipping incBSTBatch (only valid for perBatch mode)")
    else:
        if any(algo['name'] == "incNormalBST" for algo in algorithms_to_run):
            print("[INFO] Skipping incNormalBST (perEdge only)")
        algorithms_to_run = [algo for algo in algorithms_to_run if algo['name'] != "incNormalBST"]

    print(f"--- Target Dataset: {args.dataset} ---")

    # --- PHASE 0: Clean previous builds ---
    print("\n[PHASE 0: Cleaning up previous builds]")
    if not run_make_clean(MAKE_DIR):
        print("\nAborting pipeline due to an error during 'make clean'.", file=sys.stderr)
        sys.exit(1)

    # --- PHASE 1: Build the list of jobs ---
    print("\n[PHASE 1: Generating job list]")
    jobs_to_run = []
    
    if args.dataset == 'Ecoli':
        ### MODIFIED ###
        # Add a warning if -c is used with Ecoli
        if args.chromosomes:
            print("[!] Warning: '--chromosomes' argument is ignored when dataset is 'Ecoli'.")
            
        # --- ECOLI DATASET ---
        input_path = os.path.join(PROJECT_ROOT, "input", "Ecoli", "Graphs", "numeric_edgelist.txt")
        if not os.path.exists(input_path):
             print(f"[!] ERROR: Ecoli input file not found at '{input_path}'", file=sys.stderr)
             sys.exit(1)
             
        for algo in algorithms_to_run:
            output_dir = os.path.join(PROJECT_ROOT, "output", "Ecoli", "timeFile", mode_dir)
            full_output_path = os.path.join(output_dir, algo['output_filename'])
            jobs_to_run.append({
                "algo_name": algo['name'],
                "target": algo['target'],
                "input_path": input_path,
                "output_dir": output_dir,
                "full_output_path": full_output_path,
                "job_name": f"{algo['name']} on Ecoli" # Added job_name
            })
            
    elif args.dataset in ('Human', 'HPRC'):
        # --- HUMAN / HPRC DATASET (one folder per chromosome) ---
        # HPRC = Minigraph-Cactus graphs prepared by Preprocessing/mcToDag.py
        # (experimental), input/HPRC/chromosome_<X>/, same layout as Human
        
        ### MODIFIED ###
        # Define the full list of chromosomes
        ALL_HUMAN_CHROMOSOMES = [f"chromosome_{i}" for i in range(1, 23)] + ["chromosome_X", "chromosome_Y", "chromosome_M"]
        
        chromosomes_to_process = []
        
        if args.chromosomes:
            # 1. User specified a list of chromosomes
            chromosomes_to_process = args.chromosomes
            print(f"   ...Preparing jobs for {len(chromosomes_to_process)} specific human chromosomes: {', '.join(chromosomes_to_process)}")
        else:
            # 2. No specific chromosomes provided, so use the full list
            chromosomes_to_process = ALL_HUMAN_CHROMOSOMES
            print(f"   ...Preparing jobs for all {len(chromosomes_to_process)} human chromosomes.")

        # Loop over the *selected* list (either the full one or the user-specified one)
        for chr_name in chromosomes_to_process:
        ### END MODIFIED ###
            
            # This path assumes your Human graphs follow the Ecoli naming convention
            # e.g., .../input/Human/Graphs/chromosome_1_numeric_edgelist.txt
            # e.g., .../input/Human/Graphs/chromosome_X_numeric_edgelist.txt
            input_path = os.path.join(PROJECT_ROOT, "input", args.dataset, f"{chr_name}", "Graphs" ,"numeric_edgelist.txt")
            
            if not os.path.exists(input_path):
                print(f"[!] Warning: Input file not found, skipping: '{input_path}'", file=sys.stderr)
                continue # Skip this chromosome if the input doesn't exist

            for algo in algorithms_to_run:
                # This builds the output path exactly as you requested:
                # .../output/Human/chromosome_1/timeFile/
                output_dir = os.path.join(PROJECT_ROOT, "output", args.dataset, chr_name, "timeFile", mode_dir)
                full_output_path = os.path.join(output_dir, algo['output_filename'])
                
                jobs_to_run.append({
                    "algo_name": algo['name'],
                    "target": algo['target'],
                    "input_path": input_path,
                    "output_dir": output_dir,
                    "full_output_path": full_output_path,
                    "job_name": f"{algo['name']} on {chr_name}" # For logging
                })
    else:
        print(f"[!] ERROR: Unknown dataset '{args.dataset}'.", file=sys.stderr)
        sys.exit(1)

    if not jobs_to_run:
        print("\n[!] No valid jobs were generated. Check input files and dataset name.", file=sys.stderr)
        sys.exit(1)
        
    print(f"   ...Generated {len(jobs_to_run)} total jobs.")

    # --- PHASE 2: Compile and Run All Jobs (NUM_RUNS times into temp dirs) ---
    print(f"\n[PHASE 2: Compiling and Running Algorithms ({NUM_RUNS} runs for averaging)]")

    # Keep track of compiled targets so we don't re-compile
    compiled_executables = {}

    # Temporary root directory that holds one sub-folder per run
    tmp_runs_dir = os.path.join(PROJECT_ROOT, "_tmp_runs")

    for run_idx in range(1, NUM_RUNS + 1):
        print(f"\n====== RUN {run_idx}/{NUM_RUNS} ======")
        for i, job in enumerate(jobs_to_run):
            job_name = job.get("job_name", job['algo_name'])
            print(f"\n--- Processing Job {i+1}/{len(jobs_to_run)}: {job_name} ---")

            target_name = job['target']

            # Step 1: Compile only once (first run)
            if target_name not in compiled_executables:
                executable_path = compile_target_with_make(MAKE_DIR, target_name)
                if not executable_path:
                    print(f"\nAborting pipeline due to compilation error for {target_name}.", file=sys.stderr)
                    sys.exit(1)
                if not os.path.exists(executable_path):
                    print(f"[!] ERROR: Makefile produced invalid path '{executable_path}'", file=sys.stderr)
                    sys.exit(1)
                compiled_executables[target_name] = executable_path
            else:
                executable_path = compiled_executables[target_name]
                print(f"-> Target '{target_name}' already compiled.")

            # Step 2: Build a per-run temp output path that mirrors the real output path
            rel_output_dir = os.path.relpath(job['output_dir'], PROJECT_ROOT)
            tmp_output_dir = os.path.join(tmp_runs_dir, f"run_{run_idx:03d}", rel_output_dir)
            tmp_output_path = os.path.join(tmp_output_dir, os.path.basename(job['full_output_path']))
            try:
                os.makedirs(tmp_output_dir, exist_ok=True)
            except OSError as e:
                print(f"[!] ERROR: Could not create temp directory '{tmp_output_dir}': {e}", file=sys.stderr)
                sys.exit(1)

            # Step 3: Run the job, writing to the temp path
            run_start = time.time()
            if not run_cpp_program(executable_path, job['input_path'], tmp_output_path, args.mode):
                print(f"\nAborting pipeline due to execution error on job {job_name}.", file=sys.stderr)
                sys.exit(1)
            run_time = time.time() - run_start

            if PIPELINE_MODE and PIPELINE_TIMING_FILE:
                with open(PIPELINE_TIMING_FILE, 'a') as tf:
                    tf.write(f"{job_name}\t{run_time}\n")

    # --- PHASE 3: Average the NUM_RUNS outputs and write final files ---
    print(f"\n[PHASE 3: Averaging {NUM_RUNS} runs and writing final output files]")

    for job in jobs_to_run:
        rel_output_dir = os.path.relpath(job['output_dir'], PROJECT_ROOT)
        filename = os.path.basename(job['full_output_path'])

        # Collect per-run files for this job
        run_files = [
            os.path.join(tmp_runs_dir, f"run_{r:03d}", rel_output_dir, filename)
            for r in range(1, NUM_RUNS + 1)
        ]

        # Read all run files
        all_lines = []
        for rf in run_files:
            if not os.path.exists(rf):
                print(f"[!] WARNING: Expected run file missing: '{rf}'. Skipping averaging for {filename}.", file=sys.stderr)
                all_lines = []
                break
            with open(rf, 'r') as fh:
                all_lines.append(fh.readlines())

        if not all_lines:
            continue

        # Per-line averaging: numeric tokens are averaged; non-numeric tokens kept from run 1
        num_lines = len(all_lines[0])
        averaged_lines = []
        for line_idx in range(num_lines):
            tokens_per_run = [all_lines[r][line_idx].rstrip('\n').split() for r in range(NUM_RUNS)]
            template_tokens = tokens_per_run[0]
            avg_tokens = []
            for tok_idx, tok in enumerate(template_tokens):
                try:
                    values = [float(tokens_per_run[r][tok_idx]) for r in range(NUM_RUNS)]
                    avg_val = sum(values) / len(values)
                    # Preserve integer formatting when all values were whole numbers
                    if all(v == int(v) for v in values):
                        avg_tokens.append(str(int(avg_val)))
                    else:
                        avg_tokens.append(f"{avg_val:.6f}")
                except (ValueError, IndexError):
                    avg_tokens.append(tok)
            averaged_lines.append(' '.join(avg_tokens) + '\n')

        # Write to the real final output path — same as before, so downstream scripts work
        try:
            os.makedirs(job['output_dir'], exist_ok=True)
        except OSError as e:
            print(f"[!] ERROR: Could not create output directory '{job['output_dir']}': {e}", file=sys.stderr)
            sys.exit(1)

        with open(job['full_output_path'], 'w') as fh:
            fh.writelines(averaged_lines)
        print(f"   ...Averaged output written to: {job['full_output_path']}")

        # Superbubbles of the final graph (<ID>_superbubbles.txt, written by the
        # algorithm next to its timing file); identical in every run, so the
        # list of run 1 is kept.
        sb_name = filename.replace("_time.txt", "_superbubbles.txt")
        sb_run1 = os.path.join(tmp_runs_dir, "run_001", rel_output_dir, sb_name)
        if sb_name != filename and os.path.exists(sb_run1):
            shutil.copyfile(sb_run1, os.path.join(job['output_dir'], sb_name))
            print(f"   ...Superbubble list written to: {os.path.join(job['output_dir'], sb_name)}")

    # Clean up temporary run directories
    shutil.rmtree(tmp_runs_dir, ignore_errors=True)
    print(f"\n   ...Temporary run files cleaned up.")

    print("\n--- Pipeline completed successfully! ---")


if __name__ == "__main__":
    main()