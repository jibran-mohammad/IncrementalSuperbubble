import re
from collections import defaultdict
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

import os 
import sys 

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
    base_dir = str(GRAPH_DIR)
    gfa = os.path.join(base_dir, f"{REF_BASE_NAME}.gfa")
    gfa_pos = os.path.join(base_dir, f"{REF_BASE_NAME}_pos.gfa")


except FileNotFoundError as e:
    print(f"❌ Critical Error: {e}", file=sys.stderr)
    print("Please ensure this script is located somewhere inside the 'IncSuperbubble' repository.", file=sys.stderr)
    sys.exit(1)


# -------- helpers --------

WALK_TOKEN_RE = re.compile(r'([<>])s(\d+)')     # for W lines / WALK lines
EDGE_RE = re.compile(r'^EDGE\s+(s\d+)([+-])->(s\d+)([+-])\s*$')


# -------- main transforms --------

def normalize_gfa_to_positive(gfa_in, gfa_out, reverse_complement=False):
    """
    Create a positive-orientation-only GFA:
      - Every appearance of sN- in L/W is replaced by a fresh (+) node ID.
      - Original S lines are preserved; cloned S lines are appended after them.
      - Output order: H, (other non S/L/W), S (original), S (clones), L, W.
    Returns: dict neg_id("sN")->new_id("sK") used for all rewrites.
    """
    # --- First pass: capture S lines and their sequences in order ---
    s_seq = {}           # "s12" -> sequence (or "*")
    s_lines_original = []  # list of raw S lines in original order
    header_lines = []    # H lines
    other_lines = []     # any non S/L/W (e.g., comments/extensions)
    l_lines_input = []
    w_lines_input = []

    with open(gfa_in, "r") as f:
        for line in f:
            if not line.strip():
                continue
            tag = line[0]
            if tag == 'H':
                header_lines.append(line)
            elif tag == 'S':
                s_lines_original.append(line)
                parts = line.rstrip("\n").split("\t")
                node_id = parts[1]                  # e.g., s123
                seq = parts[2] if len(parts) > 2 else "*"
                s_seq[node_id] = seq
            elif tag == 'L':
                l_lines_input.append(line)
            elif tag == 'W':
                w_lines_input.append(line)
            else:
                other_lines.append(line)

    # --- Prepare mapping for negative-orientation nodes that need cloning ---
    neg_map = {}   # "s12" (when seen as '-') -> "sNEW"
    next_new_id = _next_available_numeric_id(s_seq.keys())  # find next integer > max

    def map_node(node_id, orient):
        """
        Return (mapped_node_id, mapped_orient='+') ensuring positivity.
        If orient == '+': keep node_id as-is.
        If orient == '-': clone if not seen, return clone with '+'.
        """
        nonlocal next_new_id
        if orient == '+':
            return node_id, '+'
        # need a clone
        if node_id not in neg_map:
            new_node = f"s{next_new_id}"
            neg_map[node_id] = new_node
            # create cloned S sequence (optionally reverse-complement)
            seq = s_seq.get(node_id, "*")
            if reverse_complement and seq not in ("*", "") and all(c in "ACGTNacgtn" for c in seq):
                seq = _revcomp(seq)
            s_seq[new_node] = seq  # register new node sequence
            next_new_id += 1
        return neg_map[node_id], '+'

    # --- Rewrite L lines (both endpoints become positive-only) ---
    l_lines_out = []
    for line in l_lines_input:
        parts = line.rstrip("\n").split("\t")
        # L  sA  +  sB  -
        if len(parts) < 5:
            continue
        a, ao, b, bo = parts[1], parts[2], parts[3], parts[4]
        a2, ao2 = map_node(a, '+' if ao == '+' else '-')  # normalize sign to '+'/'-'
        b2, bo2 = map_node(b, '+' if bo == '+' else '-')
        parts[1], parts[2], parts[3], parts[4] = a2, '+', b2, '+'
        l_lines_out.append("\t".join(parts) + "\n")

    # --- Rewrite W lines (all tokens become >sX with positive orientation) ---
    w_lines_out = []
    for line in w_lines_input:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 7:
            continue
        walk_field = parts[-1]
        tokens = WALK_TOKEN_RE.findall(walk_field)
        new_tokens = []
        for o, n in tokens:
            node = f"s{n}"
            orient = '+' if o == '>' else '-'
            if orient == '-' and node not in s_seq:
                raise ValueError(f"⚠️ Walk references negative {node} not defined in any S-line")
            node2, _ = map_node(node, orient)
            new_tokens.append(f">{node2}")  # always positive / '>'
        parts[-1] = "".join(new_tokens)
        w_lines_out.append("\t".join(parts) + "\n")

    # --- Build cloned S lines in deterministic order (by numeric id) ---
    cloned_ids = [neg_map[k] for k in sorted(neg_map.keys(), key=_numeric_key)]
    s_lines_clones = []
    for nid in cloned_ids:
        seq = s_seq.get(nid, "*")
        # S <id> <seq>
        s_lines_clones.append(f"S\t{nid}\t{seq}\n")

    # --- Write out in GFA order: H, other, S(original), S(cloned), L, W ---
    with open(gfa_out, "w") as out:
        out.writelines(header_lines)
        out.writelines(other_lines)
        out.writelines(s_lines_original)
        out.writelines(s_lines_clones)
        out.writelines(w_lines_out)
        out.writelines(l_lines_out)
        

    return neg_map


def _numeric_key(sid: str):
    # sid like "s123" -> 123 for sorting
    try:
        return int(sid[1:])
    except:
        return 10**12

def _next_available_numeric_id(ids):
    # ids like {"s1","s2"} -> start from max+1
    m = 0
    for sid in ids:
        if sid.startswith("s"):
            try:
                m = max(m, int(sid[1:]))
            except:
                continue
    return m + 1

def _revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(comp)[::-1]


def normalize_ordered_haplotypes_to_positive(input_txt, output_txt, neg_map):
    """
    Rewrite ordered haplotypes file:
      - EDGE sA+->sB- : replace any '-' node by its clone from neg_map; force '+'.
      - WALK >s / <s  : replace <sN with >sNEW using neg_map.
    """
    def rewrite_edge(line):
        m = EDGE_RE.match(line.strip())
        if not m:
            return line
        s, so, t, to = m.groups()
        if so == '-':
            s = neg_map.get(s, s)
            so = '+'
        if to == '-':
            t = neg_map.get(t, t)
            to = '+'
        return f"EDGE {s}{so}->{t}{to}\n"

    def rewrite_walk(line):
        # line like: WALK\t>s12<s34...
        parts = line.rstrip('\n').split('\t', 1)
        if len(parts) != 2:
            return line
        prefix, walk = parts
        rebuilt = []
        for o, n in WALK_TOKEN_RE.findall(walk):
            sid = f"s{n}"
            if o == '<':
                rebuilt.append('>' + neg_map.get(sid, sid))
            else:
                rebuilt.append('>' + sid)
        return f"{prefix}\t{''.join(rebuilt)}\n"

    with open(input_txt, 'r') as fin, open(output_txt, 'w') as fout:
        for raw in fin:
            line = raw
            if raw.startswith('EDGE'):
                line = rewrite_edge(raw)
            elif raw.startswith('WALK'):
                line = rewrite_walk(raw)
            fout.write(line)


if __name__ == "__main__":
    

    print("🔍 Normalizing GFA and haplotype files to positive orientation only...")
    # 'gfa' and 'gfa_pos' below now correctly refer to the dynamic paths
    neg_map = normalize_gfa_to_positive(gfa, gfa_pos, reverse_complement=False)
    print(f"🧩 Cloned {len(neg_map)} negative nodes. Example mappings:")
    for i, old_id in enumerate(sorted(neg_map.keys(), key=_numeric_key)):
        print(f"    {old_id} → {neg_map[old_id]}")
        if i >= 4:
            break

    print("✅ Orientation normalization complete.")
    print(f"➡️  New GFA: {gfa_pos}")
    print(f"🧩 Cloned negative nodes: {len(neg_map)}")