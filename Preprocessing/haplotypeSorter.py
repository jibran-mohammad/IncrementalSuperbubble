"""
A script to reorder haplotype data from a file.

It prioritizes haplotypes that define new structures (nodes or edges)
by placing them at the beginning of the output, followed by all
path-only (walk-only) haplotypes. This is a common preprocessing
step for building pangenome graphs consistently.
"""
import os
import pathlib
import sys

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

    # 4. --- DYNAMIC PATH CONFIGURATION ---
    # Get dataset name and ref fasta from environment variables
    DATASET_NAME = os.getenv('PIPELINE_DATASET_NAME')
    REF_FASTA_NAME = os.getenv('PIPELINE_REF_FASTA') # Not used here, but good practice

    if not DATASET_NAME:
        print("❌ Critical Error: 'PIPELINE_DATASET_NAME' environment variable not set.", file=sys.stderr)
        print("Please run this script via the main 'graphRunner.py' orchestrator.", file=sys.stderr)
        sys.exit(1)
    
    print(f"ℹ️ Configuring for dataset: {DATASET_NAME}")
    
    # Define common paths
    GRAPH_DIR = REPO_ROOT_PATH / "input" / DATASET_NAME / "Graphs"

    # --- Specific paths for this script ---
    input_filepath = str(GRAPH_DIR / "unique_elements.txt")
    output_filepath = str(GRAPH_DIR / "ordered_haplotypes.txt")

except FileNotFoundError as e:
    print(f"❌ Critical Error: {e}", file=sys.stderr)
    print("Please ensure this script is located somewhere inside the 'IncSuperbubble' repository.", file=sys.stderr)
    sys.exit(1)

def reorder_haplotypes(file_content: str) -> str:
    """
    Parses and reorders haplotype blocks from a string.

    Args:
        file_content: A string containing the entire haplotype data.

    Returns:
        A string with the reordered haplotype data.
    """
    # Split the content into blocks for each haplotype.
    # The delimiter is the start of each new haplotype record.
    # We filter out any empty strings that might result from the split.
    haplotype_blocks = file_content.strip().split('# HAPLOTYPE:')
    haplotype_blocks = [block.strip() for block in haplotype_blocks if block.strip()]

    structural_haplotypes = []
    path_only_haplotypes = []

    # Categorize each haplotype block
    for block in haplotype_blocks:
        # A block is considered "structural" if it defines new graph components.
        if 'NODE' in block or 'EDGE' in block:
            structural_haplotypes.append(block)
        else:
            path_only_haplotypes.append(block)
            
    # Reconstruct the full text for each block by adding the header back
    # The .format() method is used for clean string construction.
    structural_haplotypes_full = ["# HAPLOTYPE: {}".format(b) for b in structural_haplotypes]
    path_only_haplotypes_full = ["# HAPLOTYPE: {}".format(b) for b in path_only_haplotypes]

    # Combine the two lists
    ordered_haplotypes = structural_haplotypes_full + path_only_haplotypes_full

    # Add comments to delineate the two sections in the output
    output_parts = [
        "############################################################",
        "# SECTION 1: Haplotypes with Structural Definitions (NODE/EDGE)",
        "############################################################\n"
    ]
    output_parts.extend(structural_haplotypes_full)
    
    output_parts.append("\n############################################################")
    output_parts.append("# SECTION 2: Path-Only Haplotypes (WALK only)")
    output_parts.append("############################################################\n")
    output_parts.extend(path_only_haplotypes_full)

    # Join all parts together into a final string, with two newlines
    # between each haplotype entry for readability.
    return "\n\n".join(output_parts)


def main():
    """
    Main function to read an input file, reorder its content,
    and write to an output file using predefined paths.
    """
    # --- PLEASE EDIT THESE PATHS ---
    # Define the input and output file paths here.
    # You can use full paths (e.g., "C:/Users/YourUser/Documents/data.txt")
    # or relative paths (e.g., "data/input.txt").
    #input_filepath = "/home/dell/Research/Superbubble/Experiment/IncSuperbubble/input/Ecoli/Graphs/unique_elements.txt"
    #output_filepath = "/home/dell/Research/Superbubble/Experiment/IncSuperbubble/input/Ecoli/Graphs/ordered_haplotypes.txt"
    # -----------------------------

    try:
        print(f"-> Reading data from: {input_filepath}")
        with open(input_filepath, 'r') as f_in:
            content_from_file = f_in.read()

        # Reorder the content
        ordered_content = reorder_haplotypes(content_from_file)

        # Write to the output file
        with open(output_filepath, 'w') as f_out:
            f_out.write(ordered_content)

        print(f"-> Successfully created reordered file: {output_filepath}")

    except FileNotFoundError:
        print(f"Error: Input file '{input_filepath}' not found.")
        print("Please make sure the path and filename are correct.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


if __name__ == "__main__":
    main()

