import os
import sys
import glob
from collections import deque, defaultdict
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

    # 4. --- DYNAMIC PATH CONFIGURATION ---
    # Get dataset name and ref fasta from environment variables
    DATASET_NAME = os.getenv('PIPELINE_DATASET_NAME')
    REF_FASTA_NAME = os.getenv('PIPELINE_REF_FASTA')

    if not DATASET_NAME or not REF_FASTA_NAME:
        print("❌ Critical Error: 'PIPELINE_DATASET_NAME' or 'PIPELINE_REF_FASTA' environment variable not set.", file=sys.stderr)
        print("Please run this script via the main 'graphRunner.py' orchestrator.", file=sys.stderr)
        sys.exit(1)
    
    print(f"ℹ️ Configuring for dataset: {DATASET_NAME}")
    
    # Define common paths
    REF_BASE_NAME = os.path.splitext(REF_FASTA_NAME)[0]
    GRAPH_DIR = REPO_ROOT_PATH / "input" / DATASET_NAME / "Graphs"

    # --- Specific paths for this script ---
    file_path = str(GRAPH_DIR / f"{REF_BASE_NAME}.gfa")

except FileNotFoundError as e:
    print(f"❌ Critical Error: {e}", file=sys.stderr)
    print("Please ensure this script is located somewhere inside the 'IncSuperbubble' repository.", file=sys.stderr)
    sys.exit(1)

def parse_gfa_to_adjacency_list(file_path):
    """
    Parses a GFA file and builds a directed graph represented as an adjacency list.
    
    Each node is represented by a unique integer ID, and stores:
        - segment name
        - sequence
        - orientation (+ or -)
    
    Returns:
        tuple: (adj_list, node_info, edges)
            adj_list: dict {node_id: [neighbors]}
            node_info: dict {node_id: (segment_name, sequence, orientation)}
            edges: list of (u, v) pairs
    """
    node_info = {}
    node_map = {}     # maps (segment_name, orientation) -> node_id
    next_id = 2
    edges = []

    with open(file_path, 'r') as f:
        for line in f:
            line_stripped = line.strip()
            if not line_stripped or line_stripped.startswith('#'):
                continue
            
            parts = line_stripped.split('\t')
            line_type = parts[0]
            
            if line_type == 'S':  # Segment line
                segment_name = parts[1]
                sequence = parts[2]

                # Forward node
                fwd_key = (segment_name, '+')
                if fwd_key not in node_map:
                    node_map[fwd_key] = next_id
                    node_info[next_id] = (segment_name, sequence, '+')
                    next_id += 1

                # Reverse node
                rev_key = (segment_name, '-')
                if rev_key not in node_map:
                    node_map[rev_key] = next_id
                    node_info[next_id] = (segment_name, sequence, '-')
                    next_id += 1
            
            elif line_type == 'L':  # Link line
                if len(parts) >= 5:
                    from_seg, from_or = parts[1], parts[2]
                    to_seg, to_or = parts[3], parts[4]

                    from_key = (from_seg, from_or)
                    to_key = (to_seg, to_or)

                    if from_key not in node_map:
                        node_map[from_key] = next_id
                        node_info[next_id] = (from_seg, "N", from_or)  # "N" if no sequence given
                        next_id += 1

                    if to_key not in node_map:
                        node_map[to_key] = next_id
                        node_info[next_id] = (to_seg, "N", to_or)
                        next_id += 1

                    edges.append((node_map[from_key], node_map[to_key]))

    # Build adjacency list
    adj_list = {i: [] for i in range(2, next_id)}
    for u, v in edges:
        adj_list[u].append(v)

    return adj_list, node_info, edges


def is_dag_iterative(adj_list):
    """
    Iterative DAG check using Kahn's Algorithm.
    
    Args:
        adj_list (dict): adjacency list
    
    Returns:
        bool: True if DAG, False if cycle exists.
    """
    in_degree = defaultdict(int)
    for u in adj_list:
        for v in adj_list[u]:
            in_degree[v] += 1

    q = deque([u for u in adj_list if in_degree[u] == 0])
    visited_count = 0

    while q:
        u = q.popleft()
        visited_count += 1
        for v in adj_list[u]:
            in_degree[v] -= 1
            if in_degree[v] == 0:
                q.append(v)

    return visited_count == len(adj_list)


def check_sanitized_graphs_for_dag(cycle_dir):
    """
    Finds and checks all sanitized GFA files for DAG property.
    """
    cycle_files = sorted(glob.glob(os.path.join(cycle_dir, "*.gfa")))

    if not cycle_files:
        print(f"No  GFA files found in {cycle_dir}")
        return

    for file_path in cycle_files:
        print(f"\nChecking '{os.path.basename(file_path)}'...")
        adj_list, node_info, edges = parse_gfa_to_adjacency_list(file_path)

        if is_dag_iterative(adj_list):
            print(" -> The graph is a DAG. ✅")
        else:
            print(" -> The graph contains a cycle! ❌")


if __name__ == "__main__":
    
    if not os.path.exists(file_path):
        sys.exit(f"Error: File not found at '{file_path}'")
    
    print(f"\nChecking '{os.path.basename(file_path)}'...")
    adj_list, node_info, edges = parse_gfa_to_adjacency_list(file_path)

    if is_dag_iterative(adj_list):
        print(" -> The graph is a DAG. ✅")
        sys.exit(0)
    else:
        print(" -> The graph contains a cycle! ❌")
        sys.exit(2)