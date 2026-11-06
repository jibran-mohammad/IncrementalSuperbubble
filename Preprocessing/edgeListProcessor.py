"""
edgeListProcessor.py

Step 7 of the graph pipeline: full_graph_edgelist.txt ->
numeric_edgelist.txt (segment sX numbered X+1), the input of the
superbubble algorithms.

Not run directly: called by scripts/graphRunner.py (and, without
graphgeneration.py, by scripts/pipelineRunner.py), which pass the dataset
and reference through the environment variables PIPELINE_DATASET_NAME
(e.g. Ecoli, Human/chromosome_11) and PIPELINE_REF_FASTA. Files are in
input/<dataset>/Graphs/. See the guide at the top of scripts/graphRunner.py.

A script to process a pangenome edgelist file, converting string node IDs
to numeric IDs and adding edge counts for haplotype additions.

This script performs a two-pass process:
1.  First Pass: It reads the entire input file to discover all unique nodes
    (e.g., 's1', 's2', 's12') and creates a definitive mapping from the
    string ID to a new integer ID.
2.  Second Pass: It reads the input file again, transforming the structure
    according to the generated map and formatting rules, then writes the
    result to an output file.

Node Numbering Convention:
- A node 'sX' is converted to the integer 'X + 1'.
- This is done so that no node is represented by the number 1.
- For example: s1 -> 2, s2 -> 3, s10 -> 11.
"""
import re
import pathlib
import sys
import os

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
    input_filepath = str(GRAPH_DIR / "full_graph_edgelist.txt")
    output_filepath = str(GRAPH_DIR / "numeric_edgelist.txt")

except FileNotFoundError as e:
    print(f"❌ Critical Error: {e}", file=sys.stderr)
    print("Please ensure this script is located somewhere inside the 'IncSuperbubble' repository.", file=sys.stderr)
    sys.exit(1)

def create_node_map(filepath: str) -> dict:
    """
    First pass: Reads the edgelist file to find all unique nodes and creates a
    mapping from string ID (e.g., 's1') to a new integer ID (e.g., 2).

    Args:
        filepath: Path to the input edgelist file.

    Returns:
        A dictionary mapping each string node ID to its new integer ID.
    """
    unique_nodes = set()
    try:
        with open(filepath, 'r') as f:
            for line in f:
                # Match lines that look like edges, e.g., "s1 s2"
                # This avoids headers, comments, and count lines.
                if re.match(r'^\s*s\d+\s+s\d+\s*$', line):
                    source, target = line.strip().split()
                    unique_nodes.add(source)
                    unique_nodes.add(target)
    except FileNotFoundError:
        print(f"Error during node mapping: Input file '{filepath}' not found.")
        return {}

    # Sort nodes naturally based on the number in their name ('s1', 's2', ..., 's10')
    sorted_nodes = sorted(list(unique_nodes), key=lambda node: int(node[1:]))

    # Create the map. A node 'sX' becomes the integer 'X + 1'.
    node_map = {node_str: int(node_str[1:]) + 1 for node_str in sorted_nodes}
    return node_map

def process_and_write(input_path: str, output_path: str, node_map: dict):
    """
    Second pass: Reads the file again, processes it section by section,
    and writes the transformed, numeric output to a new file.

    Args:
        input_path: Path to the source edgelist file.
        output_path: Path to write the numeric edgelist file.
        node_map: The dictionary for converting node IDs.
    """
    with open(input_path, 'r') as f_in, open(output_path, 'w') as f_out:
        # Write a comment explaining the node numbering convention
        f_out.write("# Node Numbering Convention: Node 'sX' is represented as integer 'X+1'.\n")
        f_out.write("# e.g., s1 -> 2, s2 -> 3, s10 -> 11.\n\n")

        # The file content is processed in chunks separated by blank lines
        content_blocks = f_in.read().strip().split('\n\n')

        # --- Process Base Graph (first block) ---
        base_graph_lines = content_blocks[0].strip().splitlines()
        f_out.write(f"{base_graph_lines[0]}\n")  # Write graph name
        f_out.write(f"{base_graph_lines[1]}\n")  # Write counts
        for edge_line in base_graph_lines[2:]:
            source, target = edge_line.strip().split()
            f_out.write(f"{node_map[source]} {node_map[target]}\n")
        
        # --- Process Haplotype Blocks (remaining blocks) ---
        for hap_block in content_blocks[1:]:
            hap_lines = hap_block.strip().splitlines()
            hap_name = hap_lines[0]
            edge_lines = hap_lines[1:]

            f_out.write(f"\n{hap_name}\n")
            
            if edge_lines:
                # If the haplotype adds edges, write the count first
                f_out.write(f"{len(edge_lines)}\n")
                for edge_line in edge_lines:
                    source, target = edge_line.strip().split()
                    f_out.write(f"{node_map[source]} {node_map[target]}\n")

def main():
    """Main execution function."""
    try:
        print(f"Pass 1: Creating node map from '{input_filepath}'...")
        node_map = create_node_map(input_filepath)

        if not node_map:
            print("Node mapping could not be created. Aborting.")
            return

        print(f"Pass 2: Processing and writing to '{output_filepath}'...")
        process_and_write(input_filepath, output_filepath, node_map)
        
        print("\nProcess completed successfully.")
        print(f"Output file '{output_filepath}' has been created.")

    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}")

if __name__ == "__main__":
    main()
