import os
import subprocess
import sys
import shutil
import argparse

# --- HOW TO RUN ---
#
# This script uses command-line arguments to select the dataset and algorithms.
#
# 1. Run on Ecoli dataset (default) with ALL algorithms:
#    python3 superbubbleRunner.py
#    OR
#    python3 superbubbleRunner.py -d Ecoli
#
# 2. Run on Human dataset with ALL algorithms:
#    python3 superbubbleRunner.py -d Human
#
# 3. Run on Ecoli dataset with a SINGLE specific algorithm:
#    python3 superbubbleRunner.py -a SimpleSuperbubble
#
# 4. Run on Ecoli dataset with MULTIPLE specific algorithms:
#    python3 superbubbleRunner.py -a SimpleSuperbubble AnotherAlgorithm
#
# 5. Run on Human dataset with a SINGLE specific algorithm:
#    python3 superbubbleRunner.py -d Human -a IncrementalAlgorithm
#
# 6. Run on Human dataset with MULTIPLE specific algorithms:
#    python3 superbubbleRunner.py -d Human -a IncrementalAlgorithm AnotherAlgorithm
#
# 
# 7. Run on Human dataset for ONLY chromosome_1:
#    python3 superbubbleRunner.py -d Human -c chromosome_1
#
# 8. Run on Human dataset for chromosome_1 and chromosome_X with a specific algorithm:
#    python3 superbubbleRunner.py -d Human -a SimpleSuperbubble -c chromosome_1 chromosome_X
# 
# 9. For Batch processing: python3 superbubbleRunner.py -d Human -m perBatch

# 10.For Every-Edge processing: python3 superbubbleRunner.py -d Human -m perEdge
  
# 11.  python3 superbubbleRunner.py -d Ecoli -a SimpleSuperbubbleH1 -m perEdge 
#
# --- END HOW TO RUN ---


# --- DYNAMIC PATH CONFIGURATION ---
# Get the absolute path to this script's directory (e.g., .../IncSuperbubble/scripts)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Get the parent directory (the true project root, e.g., .../IncSuperbubble)
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

# 1. The directory where your Makefile is located.
#    The 'make' command will be run here.
MAKE_DIR = PROJECT_ROOT

PIPELINE_MODE = os.environ.get("SUPERBUBBLE_PIPELINE_MODE") == "1"

print("PIPELINE MODE:", PIPELINE_MODE)

# --- CONFIGURATION ---

# 1. A list describing each of your C++ algorithms.
#    - 'name': A friendly name for logging.
#    - 'target': The TARGET name for the Makefile (e.g., 'simple_superbubble_app').
#    - 'output_filename': The final name for the output file.
#
#    NOTE: The 'output_dir' has been removed from this list because the
#    main() function now generates it dynamically based on the dataset.
ALGORITHMS = [
    {
        "name": "SimpleSuperbubble",
        "target": "simpleSuperbubble",
        "output_filename": "SuperbubbleSimpleTimeOutput.txt"
    },

    {
        "name": "incsuperbubbleDifferentMaxO",
        "target": "incsuperbubbleDifferentMaxO",
        "output_filename": "SuperbubbleDifferentMaxO.txt"
    },

   
    #{
    #    "name": "IncrementalAlgorithm",
    #   "target": "incsuperbubble",
    #   "output_filename": "SuperbubbleIncrementalDFSOutput.txt"
    #},
    {
        "name": "AnotherAlgorithm",
        "target": "newincrementalSuperbubble",
        "output_filename": "SuperbubbleIncrementalAVLDFSOutput.txt"
    },
     
    #{
    #    "name": "SimpleSuperbubbleH1",
    #    "target": "simpleSuperbubbleh1",
    #    "output_filename": "SuperbubbleSimpleTimeOutputH1.txt"
    #},

    {
        "name": "incsuperbubblessm",
        "target": "incsuperbubblessm",
        "output_filename": "SuperbubbleIncrementalDFSOutputSSM.txt"
    },

    {
        "name": "IncrementalSuperbubbleBatch",
        "target": "incsuperbubblebatch",
        "output_filename": "SuperbubbleIncrementalBatchDFSOutput.txt"
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
                cwd=PROJECT_ROOT
            )
        print(f"   ...Finished successfully.")
        return True
    except subprocess.CalledProcessError as e:
        print(f"   [!] Execution FAILED for {os.path.basename(executable_path)}.", file=sys.stderr)
        print("--- Program error output ---", file=sys.stderr)
        print(e.stderr, file=sys.stderr)
        print("----------------------------", file=sys.stderr)
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
        choices=['Ecoli', 'Human'], # Added choices for clarity
        help="Specify which dataset to run (e.g., 'Ecoli', 'Human'). Default: Ecoli"
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
    
    args = parser.parse_args()

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
        removed = [algo['name'] for algo in algorithms_to_run if algo['name'] == "IncrementalSuperbubbleBatch"]
        algorithms_to_run = [
            algo for algo in algorithms_to_run
            if algo['name'] != "IncrementalSuperbubbleBatch"
        ]
        if removed:
            print("[INFO] Skipping IncrementalSuperbubbleBatch (only valid for perBatch mode)")

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
            
    elif args.dataset == 'Human':
        # --- HUMAN DATASET ---
        
        ### MODIFIED ###
        # Define the full list of chromosomes
        ALL_HUMAN_CHROMOSOMES = [f"chromosome_{i}" for i in range(1, 23)] + ["chromosome_X", "chromosome_y", "chromosome_M"]
        
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
            input_path = os.path.join(PROJECT_ROOT, "input", "Human", f"{chr_name}", "Graphs" ,"numeric_edgelist.txt")
            
            if not os.path.exists(input_path):
                print(f"[!] Warning: Input file not found, skipping: '{input_path}'", file=sys.stderr)
                continue # Skip this chromosome if the input doesn't exist

            for algo in algorithms_to_run:
                # This builds the output path exactly as you requested:
                # .../output/Human/chromosome_1/timeFile/
                output_dir = os.path.join(PROJECT_ROOT, "output", "Human", chr_name, "timeFile", mode_dir)
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

    # --- PHASE 2: Compile and Run All Jobs ---
    print("\n[PHASE 2: Compiling and Running Algorithms]")
    
    # Keep track of compiled targets so we don't re-compile
    compiled_executables = {}

    for i, job in enumerate(jobs_to_run):
        job_name = job.get("job_name", job['algo_name'])
        print(f"\n--- Processing Job {i+1}/{len(jobs_to_run)}: {job_name} ---")
        
        target_name = job['target']
        
        # Step 1: Compile (if we haven't already)
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

        # Step 2: Create output directory
        try:
            os.makedirs(job['output_dir'], exist_ok=True)
        except OSError as e:
            print(f"[!] ERROR: Could not create output directory '{job['output_dir']}': {e}", file=sys.stderr)
            sys.exit(1)

        # Step 3: Run the job
        if not run_cpp_program(executable_path, job['input_path'], job['full_output_path'], args.mode):
            print(f"\nAborting pipeline due to execution error on job {job_name}.", file=sys.stderr)
            sys.exit(1)
            
    print("\n--- Pipeline completed successfully! ---")


if __name__ == "__main__":
    main()