"""
haplotypeRemoval.py

Step 4 of the graph pipeline: from <reference>_pos.gfa writes the base
graph with only the reference (reduced_graph.gfa) and, per haplotype, the
nodes and edges it adds (unique_elements.txt), plus a report with the
counts (haplotype_removal_report.txt). See GRAPH_PREPARATION.md.

Rules (the graph is defined by the walks, one per genome):
  - a segment or link that no walk uses is dropped (graphgeneration.py
    already drops such segments; links are dropped here)
  - base graph = exactly the segments and links of the reference walk
  - every other segment/link is given to one haplotype: the last one, in
    W-line order, whose walk uses it (its use count reaches 0 there)
Checks (the script stops if one fails): one walk per genome, all links
+ -> + and all walk steps forward (after normalizeOrientations.py), and
the written base graph equals the reference walk with every link's
segments defined.

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
    report_file = str(GRAPH_DIR / "haplotype_removal_report.txt")
    

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
    walk_lines = 0  # number of W lines (checked against distinct sample names)

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
                walk_lines += 1
                parts = line.strip().split('\t')
                hap_name = parts[1]
                path_str = parts[-1]
                raw_walks[hap_name] = path_str
                tokens = re.findall(r'([<>])s(\d+)', path_str)  # e.g. [('>', '123'), ('<', '124')]
                # store as list of (node_id, orient_plus_minus), e.g. ('s123','+')
                walks[hap_name] = [(f"s{n}", '+' if o == '>' else '-') for (o, n) in tokens]
    return nodes, edges, walks, raw_walks, walk_lines


def check_input(edges, walks, walk_lines):
    """
    Assumptions of the rules below; stop with a clear message if one does
    not hold (e.g. for a future graph that keeps orientations).
    """
    problems = []

    # One walk per genome: walks are stored by sample name, so a second walk
    # of the same genome would silently replace the first.
    if walk_lines != len(walks):
        problems.append(f"{walk_lines} W lines but {len(walks)} distinct sample names "
                        f"(a genome with more than one walk)")

    # Links and walk steps must all be forward, as written by
    # normalizeOrientations.py; links are matched in their stored direction.
    backward_links = sum(1 for (u, uo, v, vo) in edges if (uo, vo) != ('+', '+'))
    backward_steps = sum(1 for path in walks.values() for _, o in path if o != '+')
    if backward_links or backward_steps:
        problems.append(f"{backward_links} links and {backward_steps} walk steps are not + "
                        f"(run normalizeOrientations.py first)")

    if problems:
        raise RuntimeError("haplotypeRemoval.py input check failed: " + "; ".join(problems))


def drop_unused(nodes, edges):
    """
    Drop segments and links that no walk uses (use count 0 after
    increment_counts). Without this they would never be given to a
    haplotype and would end up in the base graph, although they belong to
    no genome; a link to a segment that is not in the base graph would
    then dangle there.
    """
    unused_nodes = [n for n, d in nodes.items() if d['count'] == 0]
    unused_edges = [e for e, d in edges.items() if d['count'] == 0]

    for n in unused_nodes:
        del nodes[n]
    for e in unused_edges:
        del edges[e]

    return len(unused_nodes), len(unused_edges)


def check_base_graph(reduced_gfa, reference_path):
    """
    The written base graph must be exactly the reference walk, and every
    link in it must connect two segments it defines.
    """
    segs, links = set(), set()
    with open(reduced_gfa) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if parts[0] == 'S':
                segs.add(parts[1])
            elif parts[0] == 'L':
                links.add((parts[1], parts[2], parts[3], parts[4]))

    ref_segs = {u for u, _ in reference_path}
    ref_links = {(u, uo, v, vo) for (u, uo), (v, vo) in zip(reference_path, reference_path[1:])}

    problems = []
    dangling = sum(1 for (u, _, v, _) in links if u not in segs or v not in segs)
    if dangling:
        problems.append(f"{dangling} links to segments not in the base graph")
    if segs != ref_segs:
        problems.append(f"segments: {len(segs - ref_segs)} not on the reference walk, "
                        f"{len(ref_segs - segs)} of the reference walk missing")
    if links != ref_links:
        problems.append(f"links: {len(links - ref_links)} not on the reference walk, "
                        f"{len(ref_links - links)} of the reference walk missing")

    if problems:
        raise RuntimeError("base graph check failed: " + "; ".join(problems))

    return len(segs), len(links)


def check_acyclic(nodes, edges):
    """
    The superbubble algorithms need a DAG. Check the graph that is actually
    used (after dropping unused segments/links) with Kahn's algorithm; the
    DAG check in DAGCheck.py runs on the graph before normalization.
    """
    indegree = {n: 0 for n in nodes}
    out = defaultdict(list)
    for (u, _, v, _) in edges:
        out[u].append(v)
        indegree[v] = indegree.get(v, 0) + 1
        indegree.setdefault(u, 0)

    queue = [n for n, d in indegree.items() if d == 0]
    seen = 0
    while queue:
        n = queue.pop()
        seen += 1
        for v in out[n]:
            indegree[v] -= 1
            if indegree[v] == 0:
                queue.append(v)

    if seen != len(indegree):
        raise RuntimeError(f"graph check failed: {len(indegree) - seen} segments are on a cycle "
                           f"(the superbubble algorithms need a DAG)")


def write_report(report_file, stats):
    """Counts of this step, for checking and for describing the data."""
    with open(report_file, 'w') as f:
        for key, value in stats:
            f.write(f"{key}: {value}\n")


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
    nodes, edges, walks, raw_walks, walk_lines = parse_gfa(gfa_file)

    print(f"Nodes: {len(nodes)} | Edges: {len(edges)} | Haplotypes: {len(walks)}")

    check_input(edges, walks, walk_lines)

    stats = [("graph", gfa_file), ("walks (genomes, incl. reference)", len(walks)),
             ("reference walk", ref_sample),
             ("segments in graph", len(nodes)), ("links in graph", len(edges))]

    print("🔢 Incrementing counts...")
    increment_counts(nodes, edges, walks)

    # Segments/links on no walk belong to no genome: drop them (see drop_unused).
    unused_nodes, unused_edges = drop_unused(nodes, edges)
    print(f"🗑️  Dropped {unused_nodes} segments and {unused_edges} links that no walk uses")
    stats += [("segments used by no walk (dropped)", unused_nodes),
              ("links used by no walk (dropped)", unused_edges)]

    check_acyclic(nodes, edges)
    print(f"✅ Graph used by the algorithms is acyclic: {len(nodes)} segments, {len(edges)} links")
    stats += [("segments used (graph of all walks)", len(nodes)),
              ("links used (graph of all walks)", len(edges)),
              ("acyclic", "yes")]

    print("🧩 Identifying haplotype-unique nodes and edges...")
    process_haplotypes(nodes, edges, walks, raw_walks, unique_file, ref_sample)
    print(f"📝 Unique elements + walks saved to {unique_file}")

    print("🧹 Writing reduced graph (without unique haplotype nodes/edges)...")
    write_reduced_gfa(nodes, edges, gfa_file, reduced_gfa, ref_sample)

    base_segments, base_links = check_base_graph(reduced_gfa, walks[ref_sample])
    print(f"✅ Base graph = reference walk: {base_segments} segments, {base_links} links")

    stats += [("base graph segments (= reference walk)", base_segments),
              ("base graph links (= reference walk)", base_links)]
    added = [0, 0]
    with open(unique_file) as f:
        for line in f:
            if line.startswith("NODE"):
                added[0] += 1
            elif line.startswith("EDGE"):
                added[1] += 1
    stats += [("segments added by haplotypes", added[0]), ("links added by haplotypes", added[1])]
    write_report(report_file, stats)
    print(f"📝 Report written to {report_file}")


if __name__ == "__main__":
    main()
