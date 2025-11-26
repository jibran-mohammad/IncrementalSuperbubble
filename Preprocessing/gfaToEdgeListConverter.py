"""
A script to convert a base pangenome graph and a series of haplotype
augmentations into a consolidated edge list format.

This script performs the following steps:
1.  Parses a base graph file (in GFA-like format) to establish an initial
    set of nodes and edges.
2.  Reads a sorted haplotype file, where haplotypes defining new structures
    (NODE/EDGE) appear first.
3.  For each haplotype, it identifies any new edges that are not already present
    in the master graph. An edge can be defined explicitly with an 'EDGE' line
    or implicitly by a 'WALK' path (e.g., >s1>s2 implies edge s1->s2).
4.  Writes the entire graph representation to a single output file in a
    specific edge list format, showing the base graph first, followed by
    the incremental additions from each haplotype.
"""
import re
from typing import Set, List, Tuple, Dict, Any
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
    base_graph_filepath = str(GRAPH_DIR / "reduced_graph.gfa")
    sorted_haplotypes_filepath = str(GRAPH_DIR / "ordered_haplotypes.txt")
    output_edgelist_filepath = str(GRAPH_DIR / "full_graph_edgelist.txt")

except FileNotFoundError as e:
    print(f"❌ Critical Error: {e}", file=sys.stderr)
    print("Please ensure this script is located somewhere inside the 'IncSuperbubble' repository.", file=sys.stderr)
    sys.exit(1)

EDGE_RE = re.compile(r'^EDGE\s+(s\d+)[+-]?->(s\d+)[+-]?\s*$')

def _base(s: str) -> str:
    # keep "s123", strip a single trailing '+' or '-'
    return s[:-1] if s and s[-1] in '+-' else s

def parse_base_graph(filepath: str) -> Tuple[str, Set[str], Set[Tuple[str, str]]]:
    """
    Parses the base reference graph file to extract the name, nodes, and edges.

    Args:
        filepath: The path to the base graph file.

    Returns:
        A tuple containing:
        - The name of the graph (e.g., 'CHM13#0').
        - A set of all node IDs (e.g., {'s1', 's2', ...}).
        - A set of all edge tuples (e.g., {('s1', 's2'), ...}).
    """
    graph_name = "reference_genome"
    nodes = set()
    edges = set()
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if not parts:
                continue
            record_type = parts[0]
            if record_type == 'S':
                nodes.add(parts[1])
            elif record_type == 'L':
                source, target = parts[1], parts[3]
                edges.add((source, target))
            elif record_type == 'W'and graph_name == "reference_genome":
                graph_name = parts[1]
    return graph_name, nodes, edges

def parse_haplotype_file(filepath: str) -> List[Dict[str, Any]]:
    """
    Parses a sorted haplotype file into a list of haplotype data dictionaries.

    Args:
        filepath: The path to the sorted haplotype file.

    Returns:
        A list of dictionaries, where each dictionary represents one haplotype
        and contains its name, defined nodes, defined edges, and walk path.
    """
    parsed_haplotypes = []
    with open(filepath, 'r') as f:
        content = f.read()
    
    haplotype_blocks = content.strip().split('# HAPLOTYPE:')
    haplotype_blocks = [block.strip() for block in haplotype_blocks if block.strip()]

    for block in haplotype_blocks:
        lines = block.strip().splitlines()
        name = lines[0].strip()
        
        hap_data = {
            'name': name,
            'defined_edges': set(),
            'walk_path': None
        }
        
        for line in lines[1:]:
            if line.startswith('EDGE'):
                m = EDGE_RE.match(line.strip())
                if not m:
                    continue
                u, v = m.groups()
                hap_data['defined_edges'].add((_base(u), _base(v)))
            elif line.startswith('WALK'):
                # e.g., "WALK >s1>s2>s3"
                hap_data['walk_path'] = line.split('\t', 1)[1]
                
        parsed_haplotypes.append(hap_data)
    return parsed_haplotypes

def extract_edges_from_haplotype(hap_data: Dict[str, Any]) -> Set[Tuple[str, str]]:
    """
    Extracts all explicitly defined edges from a single haplotype's data.

    Since the input file has been pre-processed to ensure all new edges
    (including novel adjacencies) are explicitly defined with 'EDGE' lines,
    we no longer need to infer edges from the WALK path.

    Args:
        hap_data: A dictionary containing the haplotype's info.

    Returns:
        A set of all explicitly defined edges for this haplotype.
    """
    # Only return the edges that were explicitly defined in the file.
    return hap_data.get('defined_edges', set())

def main():
    """
    Main execution function.
    """
    try:
        # 1. Process the base graph
        print(f"Parsing base graph from '{base_graph_filepath}'...")
        ref_name, master_nodes, master_edges = parse_base_graph(base_graph_filepath)
        
        # 2. Process the sorted haplotypes
        print(f"Parsing haplotypes from '{sorted_haplotypes_filepath}'...")
        haplotypes = parse_haplotype_file(sorted_haplotypes_filepath)

        # 3. Build the output content
        output_lines = []
        
        # Add base graph information
        output_lines.append(ref_name)
        output_lines.append(f"{len(master_edges)} {len(master_nodes)}")
        for source, target in sorted(list(master_edges)): # Sorted for consistent output
            output_lines.append(f"{source} {target}")
            
        # 4. Process each haplotype to find and append its new edges
        for hap in haplotypes:
            hap_edges = extract_edges_from_haplotype(hap)
            
            # Find the edges that are new in this haplotype
            new_edges = hap_edges - master_edges
            # Add the haplotype header and its new edges
            output_lines.append(f"\n# HAPLOTYPE: {hap['name']}")
            if new_edges:
                
                for source, target in sorted(list(new_edges)):
                    output_lines.append(f"{source} {target}")
            
            # Update the master set of edges for the next iteration
            master_edges.update(new_edges)

        # 5. Write to the output file
        print(f"Writing final edge list to '{output_edgelist_filepath}'...")
        with open(output_edgelist_filepath, 'w') as f:
            f.write("\n".join(output_lines))
            
        print("\nProcess completed successfully.")
        print(f"Output file '{output_edgelist_filepath}' has been created.")

    except FileNotFoundError as e:
        print(f"\nError: File not found - {e.filename}")
        print("Please make sure the input files exist and the paths are correct.")
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}")


if __name__ == "__main__":
    main()

