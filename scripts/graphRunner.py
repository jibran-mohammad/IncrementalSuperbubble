#!/usr/bin/env python3
import subprocess
import sys
import pathlib
import argparse  # Import argparse
import os        # Import os

# ==============================================================================
#  Dynamic Path Discovery Function
# ==============================================================================

def find_repo_root(start_path: pathlib.Path, repo_name: str) -> pathlib.Path:
    """
    Finds the absolute path to a named repository root directory.
    
    It starts from 'start_path' (usually the script's directory)
    and walks up the filesystem tree until it finds a directory
    matching 'repo_name'.

    Args:
        start_path: The pathlib.Path object to start the search from.
                    Typically the script's own parent directory.
        repo_name: The name of the repository directory to find
                   (e.g., "IncSuperbubble").

    Returns:
        A pathlib.Path object of the found repository root.

    Raises:
        FileNotFoundError: If the repository root directory is not found
                           in any parent path.
    """
    # Start from the resolved, absolute path
    current_dir = start_path.resolve()

    # Loop until we hit the filesystem root (where parent == self)
    while current_dir != current_dir.parent:
        if current_dir.name == repo_name:
            return current_dir  # Found it!
        current_dir = current_dir.parent

    # Final check in case the root directory itself is the repo_name
    if current_dir.name == repo_name:
        return current_dir

    # If we get here, we've hit the filesystem root and found nothing
    raise FileNotFoundError(
        f"Could not find repository root '{repo_name}' in any parent path of {start_path}"
    )

# ==============================================================================
#  Script Configuration
# ==============================================================================

try:
    # 1. Define the name of the repository we are looking for
    REPO_NAME = "IncSuperbubble"
    
    # 2. Get the directory containing this script
    try:
        script_dir = pathlib.Path(__file__).resolve().parent
    except NameError:
        script_dir = pathlib.Path.cwd()

    # 3. Find the repo root path object
    REPO_ROOT_PATH = find_repo_root(script_dir, REPO_NAME)
    PREPROCESS_DIR = str(REPO_ROOT_PATH / "Preprocessing")
    INPUT_DIR = REPO_ROOT_PATH / "input" # Define base input dir
    
    print(f"ℹ️ Located repository root: {REPO_ROOT_PATH}")
    print(f"ℹ️ Set PREPROCESS_DIR to: {PREPROCESS_DIR}")

except FileNotFoundError as e:
    print(f"❌ Critical Error: {e}", file=sys.stderr)
    print("Please ensure this script is located somewhere inside the 'IncSuperbubble' repository.", file=sys.stderr)
    sys.exit(1)


# Define pipeline order
pipeline = [
    "graphgeneration.py",
    "normalizeOrientations.py",
    "DAGCheck.py",
    "haplotypeRemoval.py",
    "haplotypeSorter.py",
    "gfaToEdgeListConverter.py",
    "edgeListProcessor.py"
]

def run_script(script_name, env_vars):
    """
    Runs a subscript with a specific set of environment variables.
    """
    script_path = f"{PREPROCESS_DIR}/{script_name}"
    print(f"\n🚀 Running {script_name} for {env_vars['PIPELINE_DATASET_NAME']}...")
    
    # Get the current environment and update it with our new variables
    current_env = os.environ.copy()
    current_env.update(env_vars)
    
    try:
        # Build the command
        command = []
        if script_name == "graphgeneration.py":
            # Pass the thread count to graphgeneration.py
            command = ["python3", script_path, "-t", "20"]
        else:
            command = ["python3", script_path]

        # Run the command with the new environment
        proc = subprocess.run(
            command,
            check=True,
            env=current_env  # Pass the modified environment
        )

        if proc.returncode != 0:
            print(f"❌ {script_name} failed (rc={proc.returncode}). Stopping pipeline.")
            sys.exit(proc.returncode)
        print(f"✅ {script_name} completed successfully.")
    except Exception as e:
        print(f"❌ Error running {script_name}: {e}")
        sys.exit(1)

def get_dataset_configs(args):
    """
    Generates a list of configurations to run based on user arguments.
    Each config is a tuple: (dataset_path_name, ref_fasta_filename)
    """
    configs = [] 
    
    if args.dataset == "ecoli":
        print("Targeting 'Ecoli' dataset.")
        # Config format: (Relative path from 'input', reference fasta name)
        configs.append(("Ecoli", "GCF_000005845.2_ASM584v2_genomic.fna"))
        
    elif args.dataset == "Human":
        human_base_dir = INPUT_DIR / "Human"
        
        # Determine which chromosomes to process
        target_chromosomes = []
        if "all" in args.chromosomes:
            print("Targeting all human chromosomes.")
            # Scan the directory for valid chromosome folders
            try:
                for entry in human_base_dir.iterdir():
                    if entry.is_dir() and entry.name.startswith("chromosome_"):
                        target_chromosomes.append(entry.name)
            except FileNotFoundError:
                print(f"❌ Error: Directory not found: {human_base_dir}", file=sys.stderr)
                sys.exit(1)
        else:
            # Use the user-provided list
            target_chromosomes = [f"chromosome_{c}" for c in args.chromosomes]
            print(f"Targeting specific human chromosomes: {target_chromosomes}")

        # Create configs for each target chromosome
        for chrom_name in sorted(target_chromosomes):
            chrom_dir = human_base_dir / chrom_name
            ref_fasta = f"{chrom_name}.fa"
            
            if not (chrom_dir / ref_fasta).exists():
                print(f"⚠️ Warning: Reference fasta '{ref_fasta}' not found in '{chrom_dir}'. Skipping.")
                continue
                
            # The DATASET_NAME will be the relative path from 'input'
            dataset_path = f"Human/{chrom_name}"
            configs.append((dataset_path, ref_fasta))
            
    return configs

def main():
    """
    HOW TO RUN THIS SCRIPT:
    
    This script orchestrates the entire pangenome pipeline.
    You must run this script from your terminal (bash, zsh, etc.).
    
    Replace 'python3' with 'python' if needed for your environment.
    Replace '/path/to/Preprocessing/graphRunner.py' with the actual path to this script.
    
    --------------------------------------------------------------------------
    
    OPTION 1: Run the 'Ecoli' dataset
    
    This is the original behavior. It will process the 'Ecoli' directory.
    
    COMMAND:
    python3 /path/to/Preprocessing/graphRunner.py --dataset ecoli
    
    --------------------------------------------------------------------------
    
    OPTION 2: Run the 'Human' dataset (ALL chromosomes)
    
    This will find all directories inside 'input/Human/' that start with
    'chromosome_' (e.g., chromosome_1, chromosome_2, ..., chromosome_X)
    and run the full pipeline on each one sequentially.
    
    COMMAND:
    python3 /path/to/Preprocessing/graphRunner.py --dataset Human --chromosomes all
    
    (Note: '--chromosomes all' is the default, so you can also just run:
     python3 /path/to/Preprocessing/graphRunner.py --dataset Human)
     
    --------------------------------------------------------------------------
    
    OPTION 3: Run the 'Human' dataset (SPECIFIC chromosomes)
    
    This will run the pipeline only on the specified chromosomes.
    Just list the chromosome numbers/letters after '--chromosomes'.
    
    EXAMPLE (run on 1, 2, and X):
    python3 /path/to/Preprocessing/graphRunner.py --dataset Human --chromosomes 1 2 X
    
    EXAMPLE (run on just chromosome M):
    python3 /path/to/Preprocessing/graphRunner.py --dataset Human --chromosomes M
    
    """
    parser = argparse.ArgumentParser(description="Pangenome Processing Pipeline Runner")
    parser.add_argument(
        "-d", "--dataset",
        choices=["ecoli", "Human"],
        required=True,
        help="The dataset to process."
    )
    parser.add_argument(
        "-c", "--chromosomes",
        nargs='+',
        default=["all"],
        help="For 'human' dataset: list of chromosomes (e.g., '1' '2' 'X') or 'all'. Default is 'all'."
    )
    args = parser.parse_args()

    # Get the list of tasks to run
    dataset_configs = get_dataset_configs(args)
    
    if not dataset_configs:
        print("No valid dataset configurations found. Exiting.")
        sys.exit(1)

    print(f"\n🔗 Starting Pangenome Processing Pipeline for {len(dataset_configs)} target(s)...")

    # Loop over each dataset configuration and run the full pipeline
    for dataset_name, ref_fasta in dataset_configs:
        print(f"\n========================================================")
        print(f"Processing Dataset: {dataset_name}")
        print(f"========================================================")
        
        # Define the environment variables for this run
        env_vars = {
            "PIPELINE_DATASET_NAME": dataset_name,
            "PIPELINE_REF_FASTA": ref_fasta
        }
        
        # Run the entire pipeline for this dataset
        for script in pipeline:
            run_script(script, env_vars)
    
    print("\n🎉 All pipeline stages completed successfully for all targets!")


if __name__ == "__main__":
    main()