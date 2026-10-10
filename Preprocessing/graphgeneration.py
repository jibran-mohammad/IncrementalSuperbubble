"""
graphgeneration.py

Step 1 of the graph pipeline: builds the pangenome graph from the FASTA
files of the dataset (minichain or minigraph, chosen with
PIPELINE_GRAPH_BACKEND), aligns every genome to it and adds one walk per
genome. Writes <reference>.gfa; alignments go to input/<dataset>/GAF_Files/.
Each genome's walk is its best primary alignment record; segments on no
walk are removed. See GRAPH_PREPARATION.md, section 1.
Run by graphRunner.py as: python3 -m Preprocessing.graphgeneration -t 27

Not run directly: called by scripts/graphRunner.py (and, without
graphgeneration.py, by scripts/pipelineRunner.py), which pass the dataset
and reference through the environment variables PIPELINE_DATASET_NAME
(e.g. Ecoli, Human/chromosome_11) and PIPELINE_REF_FASTA. Files are in
input/<dataset>/Graphs/. See the guide at the top of scripts/graphRunner.py.
"""

import os
import time
import sys
import multiprocessing
import subprocess
import getopt as getopt
import re
import shlex
import shutil
import pathlib


from .graphBackend.minichain import MinichainBackend
from .graphBackend.minigraph import MinigraphBackend

from .tools import get_tools
def get_graphBackend(name, threads, tools):
    if name == "minichain":
        return MinichainBackend(
            threads=threads,
            minichain_bin=tools["minichain"],
        )
    elif name == "minigraph":
        return MinigraphBackend(
            threads=threads,
            minigraph_bin=tools["minigraph"],
        )
    else:
        raise ValueError(f"Unknown graph backend: {name}")

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
    GRAPH_BACKEND = os.getenv("PIPELINE_GRAPH_BACKEND", "minichain")
    print(f"ℹ️ Graph backend selected: {GRAPH_BACKEND}")

    if not DATASET_NAME or not REF_FASTA_NAME:
        print("❌ Critical Error: 'PIPELINE_DATASET_NAME' or 'PIPELINE_REF_FASTA' environment variable not set.", file=sys.stderr)
        print("Please run this script via the main 'graphRunner.py' orchestrator.", file=sys.stderr)
        sys.exit(1)
    
    print(f"ℹ️ Configuring for dataset: {DATASET_NAME}")
    print(f"ℹ️ Using reference fasta: {REF_FASTA_NAME}")
    
    # Define common paths
    INPUT_BASE_DIR = REPO_ROOT_PATH / "input" / DATASET_NAME
    GRAPH_DIR = INPUT_BASE_DIR / "Graphs"
    GAF_DIR = INPUT_BASE_DIR / "GAF_Files"
    
    TOOLS = get_tools(REPO_ROOT_PATH)
    for name, path in TOOLS.items():
        if not pathlib.Path(path).exists():
            print(f"⚠️ Tool not found: {name} -> {path}")

    # --- Specific paths for this script ---
    """
    VG = str(REPO_ROOT_PATH / "Pangenome_Construction/vgTool/vg/bin/./vg")
    MINICHAIN = str(REPO_ROOT_PATH / "Pangenome_Construction/minichain/./minichain")
    mash_bin = str(REPO_ROOT_PATH / "Pangenome_Construction/Mash/./mash")
    """
    graph_dir = str(GRAPH_DIR)
    gaf_dir = str(GAF_DIR)
    
    # Handle different FASTA directory structures
    if DATASET_NAME == 'Ecoli':
        fa_dir = str(INPUT_BASE_DIR / "coli3682_dataset")
    elif "Human" in DATASET_NAME:
        # For human, the DATASET_NAME *is* the fasta dir (e.g., human/chromosome_1)
        fa_dir = str(INPUT_BASE_DIR)
    elif "Virus" in DATASET_NAME:
        # Both the reference and all FASTAs are here
        fa_dir = str(INPUT_BASE_DIR / "cleanGenomes")     

    # Set the ref fasta filename (the script joins it with fa_dir later)
    ref = REF_FASTA_NAME

except FileNotFoundError as e:
    print(f"❌ Critical Error: {e}", file=sys.stderr)
    print("Please ensure this script is located somewhere inside the 'IncSuperbubble' repository.", file=sys.stderr)
    sys.exit(1)

WALK_RE = re.compile(r'[<>]s\d+')              # matches >s123 or <s456 (minichain style)
WALK_TOKEN_NUMERIC_RE = re.compile(r'([<>])s(\d+)')   # for walk-edge building & haps
WALK_P_RE = re.compile(r'(s[0-9]+)([+-])')   # matches P-line tokens like s123+ from vg
WALK_GAF_RE = WALK_RE                        # GAF walk column looks like >s123<s456...

def extract_walk_from_gaf(parts):
    """
    Return the column that actually contains a walk string like >s123<s124...
    Scans from the right so we’re resilient to column count changes.
    """
    for col in reversed(parts):
        if WALK_RE.search(col):
            return col.strip()
    return None

def count_nodes_edges(gfa_path):
    nodes = 0
    edges = 0
    with open(gfa_path, 'r') as f:
        for line in f:
            if line.startswith('S'):
                nodes += 1
            elif line.startswith('L'):
                edges += 1
    return nodes, edges

def build_reference_w_line(graph_path, ref_fasta_path, gaf_dir, threads):
    """
    Return (ref_w_line, ref_seq_id) for the reference.
    - SampleId = basename of ref_fasta (no extension)
    - HapIndex = 0
    - SeqId    = path/sequence id in graph (from vg P-line) or from GAF fallback
    - Start/End = * *
    - Walk     = >s…/<s… string
    """
    ref_sample_id = os.path.splitext(os.path.basename(ref_fasta_path))[0]
    ref_seq_id = None
    ref_walk = None

    # 1) Try vg convert (first P-line)
    cmd = f'{shlex.quote(TOOLS["vg"])} convert -g -f {shlex.quote(graph_path)} -W | grep "^P" | head -n 1'
    try:
        out = subprocess.check_output(cmd, shell=True, text=True)
        parts = out.split('\t')
        if len(parts) >= 3:
            ref_seq_id = parts[1]
            # parts[2] like: s123+,s124-,...
            toks = WALK_P_RE.findall(parts[2])
            ref_walk = ''.join(('>' if sign == '+' else '<') + seg for seg, sign in toks)
    except subprocess.CalledProcessError:
        pass

    # 2) Fallback: align the reference to graph and extract walk from GAF
    if not ref_walk:
        tmp_gaf = os.path.join(
            gaf_dir, os.path.splitext(os.path.basename(ref_fasta_path))[0] + ".ref.gaf"
        )
        align_cmd = f"{shlex.quote(TOOLS['minichain'])} -t{threads} -cx asm --vc {shlex.quote(graph_path)} {shlex.quote(ref_fasta_path)} > {shlex.quote(tmp_gaf)}"
        os.system(align_cmd)

        best = None  # (is_primary, mapq, parts)
        def _mapq(p):
            try: return int(p[11])
            except: return -1
        def _is_primary(p):
            return any(t.startswith("tp:A:P") for t in p[12:])

        with open(tmp_gaf, 'r') as gf:
            for line in gf:
                p = line.rstrip('\n').split('\t')
                if len(p) < 12:
                    continue
                cand = (_is_primary(p), _mapq(p), p)
                if best is None or (cand[0] and not best[0]) or (cand[0]==best[0] and cand[1] > best[1]):
                    best = cand

        if best:
            _, _, p = best
            ref_seq_id = p[0]  # query name (e.g., NC_004431.1)
            ref_walk = extract_walk_from_gaf(p)

    if not ref_seq_id or not ref_walk:
        raise RuntimeError("Could not derive reference W-line (missing seq id or walk).")

    ref_w_line = f"W\t{ref_sample_id}\t0\t{ref_seq_id}\t*\t*\t{ref_walk}\n"
    return ref_w_line, ref_seq_id


threads = 20


## Pass arguments
argv = sys.argv[1:]

if(len(argv)==0):
    print("help: Gen_Graph.py -h")
    sys.exit(2)

try:
    opts, args = getopt.getopt(argv, 't:')
except:
    print("usage: Gen_Graph.py -t <threads>")
    sys.exit(2)

for opt, arg in opts:
    if opt == '-h':
        print("usage: Gen_Graph.py -t <threads>")
        sys.exit()
    elif opt in ("-t"):
        threads = int(arg)
    
print("Threads : " + str(threads))

AlignPool = 4
Align_Threads = 4

if not os.path.exists(graph_dir):
    os.makedirs(graph_dir)

if not os.path.exists(gaf_dir):                                                 
    os.makedirs(gaf_dir)                                                    

# read fasta files from "Genomes" folder
fasta_files = []
for filename in os.listdir(fa_dir):
    if filename.endswith((".fa", ".fna", ".fasta")):
        if filename == ref:  
            continue    
        fasta_files.append(os.path.join(fa_dir, filename))

# `ref` is the ref filename string (e.g. "GCF_... .fna")
ref_filename = ref
ref_base = os.path.splitext(ref_filename)[0]           
ref = os.path.join(fa_dir, ref_filename)              
fasta_files.append(ref)
# add timer for execution
fasta_gfa_2 = []
data_idx = []
fasta_gfa_2.append((ref_filename, set(fasta_files)))
data_idx.append(0)
print("Number of fasta files for " + ref_filename + " : " + str(len(fasta_files)))


start_time_ = time.time()
mash_file = []
h = []
haplotypes = []
fasta_file = []
tuple_list = []
graph = []
for i in range(len(fasta_gfa_2)):
    mash_file.append(" ")
    h.append(" ")
    graph.append(" ")
    haplotypes.append([])
    fasta_file.append([])
    tuple_list.append([])


def gen_mg_gfa(idx):
    # generate graph in parallel
    h[idx] = fasta_gfa_2[idx][0]
    haplotypes[idx] = fasta_gfa_2[idx][1]

    # add timer for execution
    start_time = time.time()
    
    fasta_file[idx] = []
    for hap in haplotypes[idx]: # reference is already added
        fasta_file[idx].append(hap)

    # sort fasta_file based on mash distance
    mash_file[idx] = "mash_"+ os.path.splitext(h[idx])[0] + ".txt"
    if os.path.exists(mash_file[idx]):
        os.remove(mash_file[idx])



    # print computating mash distance
    print("Computing Mash Distance ...")
    for fasta in fasta_file[idx]:
        cmd = f"{shlex.quote(TOOLS['mash'])} dist {shlex.quote(ref)} {shlex.quote(fasta)} | grep -v '^#' | head -n 1 | cut -f 3"
        out = subprocess.check_output(cmd, shell=True)
        out = out.decode("utf-8")
        # write to file
        with open(mash_file[idx], 'a') as f:
            f.write(fasta + "\t" + out)

    fasta_file[idx] = []
    tuple_list[idx] = []
    with open(mash_file[idx], 'r') as f:
        for line in f:
            field = line.split('\t')
            fasta = field[0]
            distance = float(field[1].split('\n')[0])
            tuple_list[idx].append((distance, fasta))

    # sort tuple_list by score
    tuple_list[idx].sort(key=lambda tup : tup[0])

    # add fasta_file in sorted order
    print("Order of fasta files based on mash distance")
    for tup in tuple_list[idx]:
        print(tup[1] + "\t" + str(tup[0]))
        fasta_file[idx].append(tup[1])

    # remove mash.txt
    os.remove(mash_file[idx])
    
    graph[idx] = os.path.join(graph_dir, ref_base + ".gfa")
    
    backend = get_graphBackend(
        GRAPH_BACKEND,
        threads,
        tools=TOOLS
        )

    backend.build_graph(
        ref_fasta=ref,
        fasta_list=fasta_file[idx],
        output_gfa=graph[idx],
        label=os.path.splitext(h[idx])[0]
        )

    nodes, edges = count_nodes_edges(graph[idx])
    print(f"\n📊 Graph Statistics (raw minigraph output)")
    print(f"Nodes (S lines): {nodes}")
    print(f"Edges (L lines): {edges}\n")
    
    """
    cmd = f"{shlex.quote(MINICHAIN)} -t{threads} -cxggs {shlex.quote(ref)}"
    for file_name in fasta_file[idx]:
        if os.path.abspath(file_name) == os.path.abspath(ref):
            continue
        cmd += " " + shlex.quote(file_name) 
    cmd += " > " + shlex.quote(graph[idx])

    print('Generating Graph from ' + str(len(fasta_file[idx])) + ' files')
    
    print(f"Running command:\n{cmd}\n")


    os.system(cmd)

    print("--- Time : %s seconds ---" % (time.time() - start_time) + " for " + os.path.splitext(h[idx])[0] + ".gfa")
    """

# generate graph in parallel
num_processes = 20  # Use the number of available CPU cores
pool = multiprocessing.Pool(processes=num_processes)
pool.map(gen_mg_gfa, data_idx)
src_gfa = os.path.join(graph_dir, ref_base + ".gfa")
dst_gfa = os.path.join(graph_dir, ref_base + "_mg.gfa")
shutil.copyfile(src_gfa, dst_gfa)




# add w lines in parallel


for h, haplotypes in fasta_gfa_2:

    # add timer for execution
    start_time = time.time()
    fasta_file = []

    for hap in haplotypes: 
        fasta_file.append(hap)

    gfa_path = os.path.join(graph_dir, os.path.splitext(h)[0] + ".gfa")

    # exctract reference walk with vg
    try:
        ref_w_line, ref_id = build_reference_w_line(gfa_path, ref, gaf_dir, threads)
    except Exception as e:
        print(f"[fatal] Reference W-line build failed: {e}")
        sys.exit(2)
 
    # perform all to all alignment
    def process_fasta(fasta):
        if os.path.abspath(fasta) == os.path.abspath(ref):
            return  # skip aligning the reference; we already handled it
        # Get the base name of the input file (e.g., C4-HG002.1)
        base_name = os.path.splitext(os.path.basename(fasta))[0]

        # Construct the full path for the output GAF file
        output_gaf = os.path.join(gaf_dir, base_name + ".gaf")

        # Build the command using an f-string and the threads variable
        cmd = f"{shlex.quote(TOOLS['minichain'])} -t{Align_Threads} -cx asm --vc {shlex.quote(gfa_path)} {shlex.quote(fasta)} > {shlex.quote(output_gaf)}"
        
        # print(cmd)
        os.system(cmd)

    
    num_processes = AlignPool  # Use the number of available CPU cores
    pool = multiprocessing.Pool(processes=num_processes)
    pool.map(process_fasta, fasta_file)

    # read all gaf files
    ref_base_name = os.path.splitext(os.path.basename(ref))[0]
    gaf_files = []
    for fasta in fasta_file:
        base_name = os.path.splitext(os.path.basename(fasta))[0]
        if base_name == ref_base_name:
            continue  # skip reference .gaf if present
        gaf_path = os.path.join(gaf_dir, base_name + ".gaf")
        if os.path.exists(gaf_path):
            gaf_files.append(gaf_path)
    
    # --- generate w lines (one per GAF file, plus 1 reference) ---
    w_lines = []
    walks = []

    # push reference W once
    print("Pushing Reference Walk ...")
    w_lines.append(ref_w_line)
    walks.append(ref_id)

    for gaf in gaf_files:
        best = None  # (is_primary, mapq, parts)
        def _mapq(p):
            try: return int(p[11])
            except: return -1
        def _is_primary(p):
            return any(t.startswith("tp:A:P") for t in p[12:])

        with open(gaf, 'r') as f:
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 12:
                    continue
                qid = parts[0]
                if qid == ref_id:
                    continue  # skip the reference itself
                cand = (_is_primary(parts), _mapq(parts), parts)
                if best is None or (cand[0] and not best[0]) or (cand[0]==best[0] and cand[1] > best[1]):
                    best = cand

        if best is None:
            continue  

        _, _, parts = best
        # walk column: scan from right to be robust
        walk_col = None
        for col in reversed(parts):
            if WALK_GAF_RE.search(col):
                walk_col = col.strip()
                break
        if not walk_col:
            continue

        sample_id = os.path.splitext(os.path.basename(gaf))[0]  # one sample per GAF
        seq_id = parts[0]  

        w_lines.append(f"W\t{sample_id}\t0\t{seq_id}\t*\t*\t{walk_col}\n")
        walks.append(seq_id)

    # sanity
    print("Number of GAF Files: " + str(len(gaf_files)))
    print("Number of W Lines: " + str(len(w_lines)))
    expected_w = 1 + len(gaf_files)  # 1 reference + 1 per GAF file
    if len(w_lines) == expected_w:
        print("Contiguous Walks Generated Successfully (1 per GAF file + reference)")
    else:
        print(f"[warn] expected {expected_w} W lines, got {len(w_lines)}")    


    # append w lines to gfa
    print("Appending W Lines to GFA ...")
    with open(gfa_path, 'a') as f:
        for line in w_lines:
            f.write(line)

    # Correct graph
    s_lines = []
    vertices = []
    edges = []
    l_lines = []
    w_lines = []
    haps = {}
    walks = []
    with open(gfa_path, 'r') as f:
        for line in f:
            if line.startswith('S'):
                field = line.split('\t')
                id = field[1]
                vertices.append(id)
                s_lines.append(line)
            elif line.startswith('L'):
                field = line.split('\t')
                id_1 = field[1]
                id_2 = field[3]
                edges.append((id_1, id_2))
                l_lines.append(line)
            elif line.startswith('W'):

                parts = line.rstrip('\n').split('\t')
                if len(parts) < 7:
                    # skip malformed W lines so walks and w_lines stay in sync
                    continue

                id = parts[1]
                walks.append(id)
                w_lines.append(line)

                # Walk string is the LAST field (safer than assuming index 6)
                walk_field = parts[-1].strip()

                # Extract tokens like >s123 or <s456
                tokens = re.findall(r'([<>])s(\d+)', walk_field)

                # Record node presence per walk (ignore orientation), dedup within this walk
                seen = set()
                for _, num in tokens:
                    seg = f"s{num}"
                    if seg not in seen:
                        haps.setdefault(seg, []).append(id)
                        seen.add(seg)
    # --- Derive oriented edges from W lines and repair missing L edges ---
    # --- Check-only: derive stats from W lines; DO NOT modify L lines ---
    link_edges = set()
    for line in l_lines:
        parts = line.strip().split('\t')
        if len(parts) >= 5:
            link_edges.add((parts[1], parts[2], parts[3], parts[4]))

    walk_edges = set()
    for walk_line in w_lines:
        parts = walk_line.rstrip('\n').split('\t')
        if len(parts) < 7:
            continue
        walk_field = parts[-1].strip()
        tokens = WALK_TOKEN_NUMERIC_RE.findall(walk_field)
        for (o1, a_num), (o2, b_num) in zip(tokens, tokens[1:]):
            a = f"s{a_num}"
            b = f"s{b_num}"
            ao = '+' if o1 == '>' else '-'
            bo = '+' if o2 == '>' else '-'
            walk_edges.add((a, ao, b, bo))

    # Just report; don't add or remove any edges
    missing_in_L = walk_edges - link_edges
    print(f"[check-only] links: {len(link_edges)}, walk-edges: {len(walk_edges)}, missing_in_L: {len(missing_in_L)}")

    # Build l_lines_ by filtering existing l_lines ONLY (no additions)
    l_lines_ = []
    for line in l_lines:
        field = line.split('\t')
        a, b = field[1], field[3]
        if a in haps and b in haps:
            l_lines_.append(line)

    for v in vertices:
        if v not in haps:
            print("Error: Vertex " + v + " not found in walks and will be removed")
    
    for walk_line in w_lines:
        parts = walk_line.rstrip('\n').split('\t')
        if len(parts) < 7:
            continue
        id = parts[1]
        walk_field = parts[-1].strip()
        tokens = re.findall(r'([<>])s(\d+)', walk_field)
        for _, num in tokens:
            seg = f"s{num}"
            if id not in haps.get(seg, []):
                print(f"Error: Walk {id} not found in haps[{seg}]")
    
    # print( "haps is correct")


    # remove redundent vertices from s_lines and l_lines
    s_lines_ = []
    l_lines_ = []
    for line in s_lines:
        field = line.split('\t')
        id = field[1]
        if id in haps:
            s_lines_.append(line)
    for line in l_lines:
        field = line.split('\t')
        id_1 = field[1]
        id_2 = field[3]
        if id_1 in haps and id_2 in haps:
            l_lines_.append(line)
        

    # read every walk and check if it is valid
    for line in s_lines_:
        field = line.split('\t')
        id = field[1]
        if id in haps:
            continue
        else:
            print("Error: Vertex " + id + " not found in walks")
            sys.exit(2)
        
    for line in l_lines_:
        field = line.split('\t')
        id_1 = field[1]
        id_2 = field[3]
        if id_1 in haps and id_2 in haps:
            continue
        else:
            print("Error: Edge " + id_1 + "->" + id_2 + " not found in walks")
    
    print("Correction of Graph Done ...")

    # delete graph and rewrite
    os.remove(gfa_path)
    with open(gfa_path, 'w') as f:
        for line in s_lines_:
            f.write(line)
        for line in w_lines:
            f.write(line) 
        for line in l_lines_:
            f.write(line)

    
    # check length of walks is same as w_lines
    if len(walks) != len(w_lines):
        print("Error: walks and w_lines length not same")
        sys.exit(2)

    print("Graph : " + gfa_path + " Number of Verices : " + str(len(haps)))

    # remove all gaf files
    for gaf in gaf_files:
        os.remove(gaf)


    print("Walks for " + gfa_path + " generated successfully")
    print("--- Time : %s seconds ---" % (time.time() - start_time))

print("--- Total Graph Generation Time : %s seconds ---" % (time.time() - start_time_))
