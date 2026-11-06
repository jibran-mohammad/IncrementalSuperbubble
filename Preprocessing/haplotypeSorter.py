"""
haplotypeSorter.py

Step 5 of the graph pipeline: reads unique_elements.txt and writes the
order in which the haplotypes are added (ordered_haplotypes.txt).

Not run directly: called by scripts/graphRunner.py (and, without
graphgeneration.py, by scripts/pipelineRunner.py), which pass the dataset
and reference through the environment variables PIPELINE_DATASET_NAME
(e.g. Ecoli, Human/chromosome_11) and PIPELINE_REF_FASTA. Files are in
input/<dataset>/Graphs/. See the guide at the top of scripts/graphRunner.py.

A script to reorder haplotype data from a file.

It sorts haplotypes into three categories based on their content:
1. Haplotypes with structural changes (Edges and Nodes).
2. Haplotypes with only Node changes (no Edges).
3. Haplotypes with no changes (Walks only).
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
    Parses and reorders haplotype blocks from a string into three categories.

    Args:
        file_content: A string containing the entire haplotype data.

    Returns:
        A string with the reordered haplotype data.
    """
    # Split the content into blocks for each haplotype.
    # The delimiter is the start of each new haplotype record.
    haplotype_blocks = file_content.strip().split('# HAPLOTYPE:')
    
    # Clean up blocks (remove empty strings from split)
    haplotype_blocks = [block.strip() for block in haplotype_blocks if block.strip()]

    structural_haplotypes = []  # Contains EDGE (and likely NODE)
    node_only_haplotypes = []   # Contains NODE but NO EDGE
    no_change_haplotypes = []   # Contains neither (WALK only)

    # Categorize each haplotype block
    for block in haplotype_blocks:
        has_edge = 'EDGE' in block
        has_node = 'NODE' in block

        if has_edge:
            # Priority 1: Haplotypes with Edges (Structural)
            structural_haplotypes.append(block)
        elif has_node:
            # Priority 2: Haplotypes with only Node changes
            node_only_haplotypes.append(block)
        else:
            # Priority 3: No structural/node changes (Walk only)
            no_change_haplotypes.append(block)
            
    # Helper to reconstruct the full text for each block
    def reconstruct(blocks):
        return ["# HAPLOTYPE: {}".format(b) for b in blocks]

    structural_full = reconstruct(structural_haplotypes)
    node_only_full = reconstruct(node_only_haplotypes)
    no_change_full = reconstruct(no_change_haplotypes)

    # Build the final output list with section headers
    output_parts = []

    if structural_full:
        output_parts.append("############################################################")
        output_parts.append("# SECTION 1: Haplotypes with Structural Definitions (EDGE/NODE)")
        output_parts.append("############################################################\n")
        output_parts.extend(structural_full)
        output_parts.append("\n")

    if node_only_full:
        output_parts.append("############################################################")
        output_parts.append("# SECTION 2: Haplotypes with Only Node Changes (NODE only)")
        output_parts.append("############################################################\n")
        output_parts.extend(node_only_full)
        output_parts.append("\n")

    if no_change_full:
        output_parts.append("############################################################")
        output_parts.append("# SECTION 3: Haplotypes with No Changes (WALK only)")
        output_parts.append("############################################################\n")
        output_parts.extend(no_change_full)

    # Join all parts together into a final string.
    # Using join with "\n\n" handles spacing between entries within lists,
    # but since we added headers and explicit newlines, we can simply join.
    return "\n\n".join(output_parts).strip()


def main():
    """
    Main function to read an input file, reorder its content,
    and write to an output file using predefined paths.
    """
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