"""
haplotypeRemoval.py

Step 4 of the graph pipeline: from <reference>_pos.gfa writes the base
graph with only the reference (reduced_graph.gfa) and, per haplotype, the
nodes and edges it adds (unique_elements.txt).

Not run directly: called by scripts/graphRunner.py (and, without
graphgeneration.py, by scripts/pipelineRunner.py), which pass the dataset
and reference through the environment variables PIPELINE_DATASET_NAME
(e.g. Ecoli, Human/chromosome_11) and PIPELINE_REF_FASTA. Files are in
input/<dataset>/Graphs/. See the guide at the top of scripts/graphRunner.py.
"""
import re
from collections import defaultdict
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
    gfa_file = str(GRAPH_DIR / f"{REF_BASE_NAME}_pos.gfa")
    unique_file = str(GRAPH_DIR / "unique_elements.txt")
    reduced_gfa = str(GRAPH_DIR / "reduced_graph.gfa")
    

except FileNotFoundError as e:
    print(f"❌ Critical Error: {e}", file=sys.stderr)
    print("Please ensure this script is located somewhere inside the 'IncSuperbubble' repository.", file=sys.stderr)
    sys.exit(1)

def detect_reference_sample(gfa_path):
    """
    Detects the reference haplotype sample_id from the first W line in the GFA.
    Your graph generator always writes the reference W first.
    """
    with open(gfa_path, 'r') as f:
        for line in f:
            if line.startswith('W'):
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 7:
                    return parts[1]  # sample_id field
    raise RuntimeError("No W line found to detect reference sample.")


def parse_gfa(file_path):
    nodes = defaultdict(lambda: {'seq': '', 'count': 0})
    edges = defaultdict(lambda: {'src_orient': '','dest_orient': '','info': '','count': 0})
    walks = defaultdict(list)
    raw_walks = {}  # To preserve original walk string

    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('S'):
                parts = line.strip().split('\t')
                node_id = parts[1]
                seq = parts[2] if len(parts) > 2 else ''
                nodes[node_id]['seq'] = seq
            elif line.startswith('L'):
                parts = line.strip().split('\t')
                src, src_orient = parts[1], parts[2]
                dest, dest_orient = parts[3], parts[4]
                edges[(src, src_orient, dest, dest_orient)] = {
                        'src_orient': src_orient,
                        'dest_orient': dest_orient,
                        'info': line.strip(),
                        'count': 0
                }
            elif line.startswith('W'):
                parts = line.strip().split('\t')
                hap_name = parts[1]
                path_str = parts[-1]
                raw_walks[hap_name] = path_str
                tokens = re.findall(r'([<>])s(\d+)', path_str)  # e.g. [('>', '123'), ('<', '124')]
                # store as list of (node_id, orient_plus_minus), e.g. ('s123','+')
                walks[hap_name] = [(f"s{n}", '+' if o == '>' else '-') for (o, n) in tokens]
    return nodes, edges, walks, raw_walks


def increment_counts(nodes, edges, walks):
    for hap, path in walks.items():
        # path: list of (node_id, orient)
        for i, (u, uo) in enumerate(path):
            nodes[u]['count'] += 1
            if i + 1 < len(path):
                v, vo = path[i+1]
                key = (u, uo, v, vo)
                if key in edges:
                    edges[key]['count'] += 1


def process_haplotypes(nodes, edges, walks, raw_walks, unique_file, ref_sample):
    with open(unique_file, 'w') as f:
        for hap, path in walks.items():
            #Skip the reference haplotype CHM13#0
            if hap == ref_sample:
                print(f"Skipping reference haplotype: {hap}")
                continue
            unique_nodes = []
            unique_edges = []

            for i in range(len(path)):
                u, uo = path[i]         # ← unpack (node_id, orient)
                if u in nodes:
                    nodes[u]['count'] -= 1
                    if nodes[u]['count'] == 0:
                        unique_nodes.append(u)

                if i < len(path) - 1:
                    (u, uo) = path[i]
                    (v, vo) = path[i+1]
                    edge_key = (u, uo, v, vo)
                    if edge_key in edges:
                        edges[edge_key]['count'] -= 1
                        if edges[edge_key]['count'] == 0:
                            unique_edges.append(edge_key)

            # Write results for this haplotype
            f.write(f"# HAPLOTYPE: {hap}\n")
            for u in unique_nodes:
                f.write(f"NODE\t{u}\n")
                del nodes[u]
            for (u, uo, v, vo) in unique_edges:
                f.write(f"EDGE\t{u}{uo}->{v}{vo}\n")
                del edges[(u, uo, v, vo)]
            f.write(f"WALK\t{raw_walks[hap]}\n\n")


def write_reduced_gfa(nodes, edges, original_gfa, output_gfa, ref_sample):
    with open(original_gfa, 'r') as infile, open(output_gfa, 'w') as outfile:
        for line in infile:
            if line.startswith('S'):
                parts = line.strip().split('\t')
                if parts[1] in nodes:
                    outfile.write(line)
            elif line.startswith('L'):
                parts = line.strip().split('\t')
                src, src_orient, dest, dest_orient = parts[1], parts[2], parts[3], parts[4]
                if (src, src_orient, dest, dest_orient) in edges:
                    outfile.write(line)
            elif line.startswith('W'):
                 # Only write the reference haplotype walk
                 if line.split('\t')[1] == ref_sample:
                    outfile.write(line)
            else:
                outfile.write(line)
    print(f"✅ Reduced GFA written to {output_gfa}")


def main():
     
    print("🔎 Detecting reference haplotype from first W line...")
    ref_sample = detect_reference_sample(gfa_file)
    print(f"Reference haplotype detected: {ref_sample}")

    print("🔍 Parsing GFA...")
    nodes, edges, walks, raw_walks = parse_gfa(gfa_file)

    print(f"Nodes: {len(nodes)} | Edges: {len(edges)} | Haplotypes: {len(walks)}")

    print("🔢 Incrementing counts...")
    increment_counts(nodes, edges, walks)

    print("🧩 Identifying haplotype-unique nodes and edges...")
    process_haplotypes(nodes, edges, walks, raw_walks, unique_file, ref_sample)
    print(f"📝 Unique elements + walks saved to {unique_file}")

    print("🧹 Writing reduced graph (without unique haplotype nodes/edges)...")
    write_reduced_gfa(nodes, edges, gfa_file, reduced_gfa, ref_sample)


if __name__ == "__main__":
    main()
