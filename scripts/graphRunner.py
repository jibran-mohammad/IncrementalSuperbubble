#!/usr/bin/env python3
import subprocess
import sys
import pathlib

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
    print(f"ℹ️ Located repository root: {REPO_ROOT_PATH}")

    # 4. --- YOUR MANUAL PATH ---
    # Manually construct the PREPROCESS_DIR using the found root path
    # You can use / from pathlib or os.path.join
    PREPROCESS_DIR = str(REPO_ROOT_PATH / "Preprocessing")
    
    # Or, if you prefer using os.path.join:
    # import os
    # PREPROCESS_DIR = os.path.join(str(REPO_ROOT_PATH), "Preprocessing")

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

def run_script(script_name):
    script_path = f"{PREPROCESS_DIR}/{script_name}"
    print(f"\n🚀 Running {script_name} ...")
    try:
        if script_name == "graphgeneration.py":
            proc = subprocess.run(
                ["python3", script_path, "-t", "20"],  # or any number of threads you want
                check=True
            )
        else:
            proc = subprocess.run(
                ["python3", script_path],
                check=True
            )

        if proc.returncode != 0:
            print(f"❌ {script_name} failed (rc={proc.returncode}). Stopping pipeline.")
            sys.exit(proc.returncode)
        print(f"✅ {script_name} completed successfully.")
    except Exception as e:
        print(f"❌ Error running {script_name}: {e}")
        sys.exit(1)

def main():
    print("🔗 Starting Pangenome Processing Pipeline...\n")
    for script in pipeline:
        run_script(script)
    print("\n🎉 All pipeline stages completed successfully!")


if __name__ == "__main__":
    main()
