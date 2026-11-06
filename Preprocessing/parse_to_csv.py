"""
parse_to_csv.py

Merges the timing files of all algorithms of one dataset and mode into one
CSV, haplotype_build_comparison.csv, which the plot and table scripts read
(plotcsv.py, csvToLatex.py, datasetWiseCaseTables.py,
plotChromosomeComparison*.py).

=====================================================================
HOW TO RUN
=====================================================================

Run from the repository root, with the virtual environment active
(source venv/bin/activate), after scripts/superbubbleRunner.py has
written the timing files (or after copying them from the server into
the same output/ folders).

    # E. coli, per-edge results
    python3 Preprocessing/parse_to_csv.py -d ecoli -m peredge

    # E. coli, per-haplotype results
    python3 Preprocessing/parse_to_csv.py -d ecoli -m perbatch

    # human chromosomes 1, 2 and 22
    python3 Preprocessing/parse_to_csv.py -d human -c 1 2 22 -m peredge

    # all human chromosomes (every output/Human/chromosome_* folder)
    python3 Preprocessing/parse_to_csv.py -d human -m peredge

    # E. coli and all human chromosomes
    python3 Preprocessing/parse_to_csv.py -d all -m peredge

Options:

    -d / --dataset      ecoli | human | all (lowercase)               (required)
    -c / --chromosome   human only: chromosome numbers, e.g. 1 2 22   (default: all)
    -m / --mode         peredge | perbatch (lowercase)                (default: perbatch)

=====================================================================
INPUT / OUTPUT
=====================================================================

Input, per folder output/Ecoli/timeFile/<perEdge|perbatch>/ or
output/Human/chromosome_<X>/timeFile/<perEdge|perbatch>/:

    static_time.txt, incRecompute_time.txt, incNormalBST_time.txt, incBST_time.txt,
    incMinIBST_time.txt, incMinIBSTH1_time.txt, incRecomputeBatch_time.txt,
    incBSTBatch_time.txt, incBSTBatchH1_time.txt, incBSTBatchH2_time.txt
    (written by superbubbleRunner.py; missing files are skipped with a
    warning, e.g. incBSTBatch_time.txt in perEdge mode)

Output, in the same folder:

    haplotype_build_comparison.csv
        one row per update (Step), with HaplotypeName, Nodes, Edges,
        SuperbubbleCount (of static), the per-algorithm columns
        <ID><Metric> (e.g. incMinIBSTTotalTime, incBSTAffectedVertices;
        IDs: see ALGORITHMS.md) and the cases Case1..Case5

Next steps: python3 Preprocessing/plotcsv.py ..., python3 Preprocessing/csvToLatex.py ...
"""

import pandas as pd
import re
import os
import pathlib
import sys
import argparse

# ==============================================================================
#  Dynamic Path Discovery Function (from graphRunner.py)
# ==============================================================================

def find_repo_root(start_path: pathlib.Path, repo_name: str) -> pathlib.Path:
    """
    Finds the absolute path to a named repository root directory.
    
    It starts from 'start_path' (usually the script's directory)
    and walks up the filesystem tree until it finds a directory
    matching 'repo_name'.
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
#  Helper Function for Time Conversion
# ==============================================================================

def robust_float_convert(val_str: str) -> float:
    """
    Converts a string to a float, correctly handling
    standard decimal or scientific notation (e.g., "1.23e-05").
    Returns 0.0 on conversion error.
    """
    try:
        return float(val_str)
    except ValueError:
        print(f"    Warning: Could not convert '{val_str}' to float. Using 0.0.")
        return 0.0
    except TypeError:
        print(f"    Warning: Invalid type for float conversion. Using 0.0.")
        return 0.0

# ==============================================================================
#  Parsing Function (Updated to handle new file format)
# ==============================================================================

def parse_file(filename: pathlib.Path, time_column_name: str, mode: str):
    data = []
    peak_memory = 0
    if not filename.exists():
        print(f"    Warning: File not found: {os.path.basename(filename)}. Skipping.")
        return pd.DataFrame(), 0, {}

    print(f"    Parsing {filename.name}...")
    with open(filename, 'r') as f:
        lines = f.readlines()

    line_iter = iter(lines)
    
    # Define metric columns
    BASE_METRIC_COLUMNS = ['UpdateDS_Time', 'DFS_Time', 'TopologicalOrderTime', 'MinI_Time', 'MaxO_Time', 'SB_Compute_Time', 'Total_Time', 'SB_Count', 'AffectedVertices', 'ExploredVertices', 'TotalAffectedEdges', 'TotalExploredEdges']
    
    # Timing files are named <ID>_time.txt (see ALGORITHMS.md); the
    # incremental algorithms write case columns, and some also explored counts.
    algorithm_id = filename.name[:-len("_time.txt")]

    is_incremental_file = algorithm_id in {"incRecompute", "incNormalBST", "incBST", "incMinIBST", "incMinIBSTH1",
                                           "incRecomputeBatch", "incBSTBatch", "incBSTBatchH1", "incBSTBatchH2"}

    has_explored = algorithm_id in {"incNormalBST", "incBST", "incMinIBST", "incMinIBSTH1",
                                    "incRecomputeBatch", "incBSTBatch", "incBSTBatchH1", "incBSTBatchH2"}
    if is_incremental_file:
        CASE_COLUMNS = ['Case1', 'Case2', 'Case3', 'Case4', 'Case5']
        if has_explored:
            EXTRA_COLUMNS = ['AffectedVertices', 'ExploredVertices', 'TotalAffectedEdges', 'TotalExploredEdges']
        else:
            EXTRA_COLUMNS = ['AffectedVertices', 'TotalAffectedEdges']

        FULL_METRIC_COLUMNS = BASE_METRIC_COLUMNS[:-4] + CASE_COLUMNS + EXTRA_COLUMNS + ['TotalVertices', 'TotalEdges']
        EXPECTED_DATA_PARTS_LEN = 2 + len(FULL_METRIC_COLUMNS)
    else:
        CASE_COLUMNS = []
        FULL_METRIC_COLUMNS = [
            'UpdateDS_Time',
            'DFS_Time',
            'TopologicalOrderTime',
            'MinI_Time',
            'MaxO_Time',
            'SB_Compute_Time',
            'Total_Time',
            'SB_Count'
        ]
        EXPECTED_DATA_PARTS_LEN = 2 + len(FULL_METRIC_COLUMNS)
        
    # Determine PREFIX and SOURCE METRIC
    if filename.name == "static_time.txt":
        source_time_metric = 'Total_Time'
        prefix = 'static_'
    elif filename.name == "incRecompute_time.txt":
        source_time_metric = 'DFS_Time'
        prefix = 'incRecompute_'
    elif filename.name == "incNormalBST_time.txt":
        source_time_metric = 'DFS_Time'
        prefix = 'incNormalBST_'
    elif filename.name == "incBST_time.txt":
        source_time_metric = 'DFS_Time'
        prefix = 'incBST_'
    elif filename.name == "incBSTBatch_time.txt":
        source_time_metric = 'DFS_Time'
        prefix = 'incBSTBatch_' 
    elif filename.name == "incRecomputeBatch_time.txt":
        source_time_metric = 'DFS_Time'
        prefix = 'incRecomputeBatch_'
    elif filename.name == "incBSTBatchH1_time.txt":
        source_time_metric = 'DFS_Time'
        prefix = 'incBSTBatchH1_' 
    elif filename.name == "incBSTBatchH2_time.txt":
        source_time_metric = 'DFS_Time'
        prefix = 'incBSTBatchH2_'
    elif filename.name == "incMinIBST_time.txt":
        source_time_metric = 'DFS_Time'   
        prefix = 'incMinIBST_'       
    elif filename.name == "incMinIBSTH1_time.txt":
        source_time_metric = 'DFS_Time'
        prefix = 'incMinIBSTH1_'

    current_haplotype_edges = 0
    current_haplotype_nodes = 0

    baseline_memory = 0
    pending_block_rows = []
    prev_edges = 0
    prev_nodes = 0
    genome_names = {}            # every genome of the file, also those that add nothing

    def flush(after_edges, after_nodes):
        """The programs write a genome's rows first and then its size line
        ("Edges Nodes" AFTER the genome), so the rows get the size that
        follows them. A genome that added nothing has no size line of its
        own; its rows are flushed with the next genome's and keep the size
        before (prev_*), since the graph did not change."""
        if not pending_block_rows:
            return
        last = pending_block_rows[-1]['Haplotype_Name']
        for row in pending_block_rows:
            own = row['Haplotype_Name'] == last
            row['Edges'] = after_edges if own else prev_edges
            row['Nodes'] = after_nodes if own else prev_nodes
            data.append(row)
        pending_block_rows.clear()

    try:
        # --- Lines 1-3: Initialization and Header Skip ---
        stats_line = next(line_iter).strip()
        if "PeakMemory" in stats_line: return pd.DataFrame(), 0, {}
        stats_parts = stats_line.split()
        if len(stats_parts) < 2: return pd.DataFrame(), 0, {}
        initial_edges = int(stats_parts[0])
        initial_nodes = int(stats_parts[1])
        current_haplotype_edges = initial_edges
        current_haplotype_nodes = initial_nodes

        build_line = next(line_iter).strip()
        if not build_line.startswith("BaseGraph_Build"): raise ValueError(f"Expected 'BaseGraph_Build' but got: {build_line}")
        build_parts = build_line.split()
        base_time = robust_float_convert(build_parts[1])
        base_superbubbles = int(build_parts[2]) # The starting total count (e.g., 5602)
        
        base_entry = {
            'Step': 0, 
            'Haplotype_Name': 'BaseGraph_Build',
            'Nodes': initial_nodes,
            'Edges': initial_edges,
            f'{prefix}Time_s': base_time 
        }
        
        # Populate all prefixed metric columns
        for col in FULL_METRIC_COLUMNS:
             if col == 'SB_Count':
                 # Map BaseGraph total count to the prefixed SB_Count column
                 base_entry[f'{prefix}{col}'] = base_superbubbles
             else:
                 base_entry[f'{prefix}{col}'] = 0.0 if 'Time' in col or col == 'SB_Compute_Time' else 0

        data.append(base_entry)

        try:
            # The column header follows the BaseGraph_Build line: skip it
            # (and blank lines), but nothing more. (It used to read one more
            # line whenever the header was not a number, i.e. always, and
            # dropped the first update: the first edge in perEdge mode, the
            # first genome in perBatch mode.)
            header_line = next(line_iter).strip()
            while not header_line:
                header_line = next(line_iter).strip()
            # Read the columns by the names in this header: the programs
            # write them in their own order (e.g. ... TotalVertices TotalEdges
            # TotalAffectedEdges TotalExploredEdges), which is not the order of
            # the fixed lists above (reading by position put the values into
            # the wrong columns). The fixed lists stay as the fallback for a
            # file without a usable header.
            header_names = header_line.split()
            # files written before the rename call the first step Insert_Time
            header_names = ['UpdateDS_Time' if h == 'Insert_Time' else h for h in header_names]
            if (len(header_names) > 2 and header_names[1] == 'Haplotype_Name'
                    and 'Total_Time' in header_names):
                FULL_METRIC_COLUMNS = header_names[2:]
                EXPECTED_DATA_PARTS_LEN = len(header_names)
            else:
                print(f"    Warning: no column header in {filename.name}; reading the columns by position")
        except StopIteration:
            pass
        
        prev_edges = initial_edges
        prev_nodes = initial_nodes
        global_step = 0 
        # --- Data Lines ---
        while True:
            # --- Line N: Edges and Nodes ---
            stats_line = next(line_iter).strip()
            
            if stats_line.startswith("BaselineMemory"):
                # the program's memory before reading the input (baseline of PeakMemory)
                baseline_memory = int(stats_line.split(':')[1].strip())
                continue
            if "PeakMemory" in stats_line:
                # Extract the peak memory value
                peak_memory = int(stats_line.split(':')[1].strip())
                # rows left without a size line (genomes that added nothing at
                # the end): the graph did not change
                flush(prev_edges, prev_nodes)
                break
            if not stats_line: continue

            stats_parts = stats_line.split()

            # Both modes: the file has the update rows of a genome first and
            # then its "Edges Nodes" line, so perBatch is read like perEdge
            # (the old perBatch branch below expected "Edges Nodes" before the
            # row and is no longer used).
            if False and mode == "perbatch":

                # stats_line MUST be "Edges Nodes"
                if not (len(stats_parts) == 2 and stats_parts[0].isdigit() and stats_parts[1].isdigit()):
                    continue  # ignore noise / wrong lines safely

                current_haplotype_edges = int(stats_parts[0])
                current_haplotype_nodes = int(stats_parts[1])

                # Now the next line is the real haplotype build line
                try:
                    build_line = next(line_iter).strip()
                except StopIteration:
                    break

                    
                
                if "PeakMemory" in build_line: 
                    peak_memory = int(build_line.split(':')[1].strip())
                    break
                if not build_line: continue

                build_parts = build_line.split()
                
                # --- Handle short zero-time lines ---
                is_short_zero_entry = (
                    len(build_parts) == 4 and 
                    build_parts[0].isdigit() and 
                    robust_float_convert(build_parts[2]) == 0.0
                )
                
                if len(build_parts) != EXPECTED_DATA_PARTS_LEN and not is_short_zero_entry:
                     continue 

                if is_short_zero_entry:
                    step = int(build_parts[0])
                    haplotype_name = build_parts[1]
                    sb_count = int(build_parts[3])
                    
                    metrics = [0] * len(FULL_METRIC_COLUMNS)
                    metrics[FULL_METRIC_COLUMNS.index('SB_Count')] = sb_count 
                
                else: 
                    step = int(build_parts[0])
                    haplotype_name = build_parts[1]
                    
                    metrics = []
                    for i, val in enumerate(build_parts[2:]):
                        col_name = FULL_METRIC_COLUMNS[i]
                        if 'Time' in col_name or col_name == 'SB_Compute_Time':
                            metrics.append(robust_float_convert(val))
                        else:
                            try:
                                metrics.append(int(val))
                            except ValueError:
                                 metrics.append(0)

                entry = {
                    'Step': step,
                    'Haplotype_Name': haplotype_name,
                    'Nodes': current_haplotype_nodes,
                    'Edges': current_haplotype_edges,
                }
                
                # Map metrics with the prefix
                for name, value in zip(FULL_METRIC_COLUMNS, metrics):
                    entry[f'{prefix}{name}'] = value

                # Assign the dynamic comparison time column
                entry[f'{prefix}Time_s'] = entry[f'{prefix}{source_time_metric}']

                data.append(entry)
            else:
                # =================== PER-EDGE MODE ===================
                # In this mode file format is:
                #
                #   Edges Nodes
                #   Edge_Number Haplotype_Name ...
                #   <edge rows...>
                #   Edges Nodes
                #   <edge rows...>
                #   ...
                #   PeakMemory ...
                #
                # Important: an "Edges Nodes" line is the size AFTER the
                # rows before it (their genome). So we keep the rows in a
                # buffer and flush them with the next "Edges Nodes" line (or
                # at PeakMemory/EOF, with the last size).

                # (1) If this is an "Edges Nodes" line
                if len(stats_parts) == 2 and stats_parts[0].isdigit() and stats_parts[1].isdigit():
                    new_edges = int(stats_parts[0])
                    new_nodes = int(stats_parts[1])

                    # Flush the previously collected edge rows (if any)
                    # the size line belongs to the rows before it (the
                    # size after their genome); it used to be given to the
                    # rows after it, one genome late
                    flush(new_edges, new_nodes)

                    # Update the "previous block size" to this new metadata
                    prev_edges = new_edges
                    prev_nodes = new_nodes

                    # Move on to next line
                    continue


                # (3) Otherwise, treat it as an edge-update row
                build_parts = stats_parts

                is_short_zero_entry = (
                    len(build_parts) == 4 and
                    build_parts[0].isdigit() and
                    robust_float_convert(build_parts[2]) == 0.0
                )

                if len(build_parts) != EXPECTED_DATA_PARTS_LEN and not is_short_zero_entry:
                    # ignore noise lines safely
                    continue

                # edge number in the file (YOU WANT THIS EXACTLY)
                local_edge_no = int(build_parts[0])
                global_step += 1
                # perEdge: running number over all edges; perBatch: the
                # update (genome) number written by the algorithm
                step = local_edge_no if mode == "perbatch" else global_step
                haplotype_name = build_parts[1]

                # Parse metrics
                if is_short_zero_entry:
                    # edgeNo hapName 0 SB_Count
                    sb_count = int(build_parts[3])

                    metrics = [0] * len(FULL_METRIC_COLUMNS)
                    metrics[FULL_METRIC_COLUMNS.index('SB_Count')] = sb_count
                    
                else:
                    metrics = []
                    for i, val in enumerate(build_parts[2:]):
                        col_name = FULL_METRIC_COLUMNS[i]
                        if 'Time' in col_name or col_name == 'SB_Compute_Time':
                            metrics.append(robust_float_convert(val))
                        else:
                            try:
                                metrics.append(int(val))
                            except ValueError:
                                metrics.append(0)

                # Buffer the row; Nodes/Edges will be filled when we see the next "Edges Nodes"
                entry = {
                    'Step': step,                 # GLOBAL step (unique in CSV)
                    'LocalEdgeNo': local_edge_no, # edge number in that haplotype block
                    'Haplotype_Name': haplotype_name,
                    'Nodes': 0,
                    'Edges': 0,
                }

                for name, value in zip(FULL_METRIC_COLUMNS, metrics):
                    entry[f'{prefix}{name}'] = value

                entry[f'{prefix}Time_s'] = entry[f'{prefix}{source_time_metric}']

                pending_block_rows.append(entry)
                genome_names[haplotype_name] = True
                    
            
    except StopIteration:
        flush(prev_edges, prev_nodes)
    except Exception as e:
        print(f"    ❌ Error parsing {filename}: {e}", file=sys.stderr)

    df = pd.DataFrame(data)
    if not df.empty:
        df['Step'] = df['Step'].astype(object)
        df.loc[df['Step'] == 0, 'Step'] = 'Base'
    
    info = {'genomes': len(genome_names), 'final_edges': prev_edges, 'final_nodes': prev_nodes,
            'baseline_memory': baseline_memory}
    return df, peak_memory, info

def real_graph_size(time_dir_path: pathlib.Path):
    """Vertices and edges of the input graph without the artificial source 0,
    sink 1 and their edges: the distinct vertices and edges of
    input/<dataset>/Graphs/numeric_edgelist.txt (the graph after the last
    genome). None if the edge list is not found."""
    parts = time_dir_path.resolve().parts
    if "output" not in parts:
        return None
    i = len(parts) - 1 - parts[::-1].index("output")
    dataset = parts[i + 1:parts.index("timeFile", i)]
    edge_list = pathlib.Path(*parts[:i]) / "input" / pathlib.Path(*dataset) / "Graphs" / "numeric_edgelist.txt"
    if not edge_list.exists():
        print(f"  Warning: {edge_list} not found; no real graph size")
        return None
    edges, header = set(), 0
    with open(edge_list) as f:
        for line in f:
            p = line.split()
            if not p or p[0].startswith("#"):
                continue
            header += 1
            if header <= 2:            # graph name, then "edges nodes" of the base graph
                continue
            if len(p) == 2:
                edges.add((p[0], p[1]))
    vertices = {v for e in edges for v in e}
    return len(vertices), len(edges)

# ==============================================================================
#  New Processing Function
# ==============================================================================

def process_time_directory(time_dir_path: pathlib.Path, mode: str):
    """
    Parses all three time files, prefixes the common metrics, and merges them.
    """
    # ADD THIS CHECK:
    required = ["static_time.txt", 
                "incRecompute_time.txt", 
                "incNormalBST_time.txt",
                "incBST_time.txt",
                "incBSTBatch_time.txt",
                "incMinIBST_time.txt",
                "incMinIBSTH1_time.txt",
                "incRecomputeBatch_time.txt",
                "incBSTBatchH1_time.txt",
                "incBSTBatchH2_time.txt"]
    
    dir_name = time_dir_path.name.lower()

    if mode == "peredge" and "peredge" not in dir_name:
        return

    if mode == "perbatch" and "perbatch" not in dir_name:
        return


    # Only proceed if AT LEAST one of these files exists in the current folder
    if not any((time_dir_path / f).exists() for f in required):
        return
    


    print(f"\nProcessing directory: {time_dir_path}")
    
    if not time_dir_path.is_dir():
        print(f"  Warning: Directory not found. Skipping.")
        return

    # Define files, corresponding time columns, and their prefixes
    ALGORITHMS = {
        "static": "static_time.txt",
        "incRecompute": "incRecompute_time.txt",
        "incNormalBST": "incNormalBST_time.txt",
        "incBST": "incBST_time.txt",
        "incBSTBatch": "incBSTBatch_time.txt",
        "incMinIBST": "incMinIBST_time.txt",
        "incMinIBSTH1": "incMinIBSTH1_time.txt",
        "incRecomputeBatch": "incRecomputeBatch_time.txt",
        "incBSTBatchH1": "incBSTBatchH1_time.txt",
        "incBSTBatchH2": "incBSTBatchH2_time.txt"
    }

    # perEdge-only algorithms (PEREDGE_ONLY in superbubbleRunner.py): in
    # perbatch the comparison is static, incMinIBSTH1 and the batch
    # algorithms; perbatch time files of the others (from earlier runs) are
    # not read
    if mode == "perbatch":
        for a in ("incNormalBST", "incRecompute", "incBST", "incMinIBST"):
            ALGORITHMS.pop(a)

    files_to_parse = {
        time_dir_path / fname: f"{prefix}_Time_s"
        for prefix, fname in ALGORITHMS.items()
    }

    all_data = []
    peak_memory_kb = {}        # algorithm ID -> PeakMemory(KB) of its time file
    baseline_memory_kb = {}    # algorithm ID -> BaselineMemory(KB) (0 if the file has none)
    dataset_info = {}          # genomes and final graph size (the same in every time file)

    for (prefix, _fname), (filepath, time_col) in zip(ALGORITHMS.items(), files_to_parse.items()):
        df, peak_mem, info = parse_file(filepath, time_col, mode)
        if info and not dataset_info:
            dataset_info = info
        if not df.empty:
            all_data.append(df)
            if peak_mem:
                peak_memory_kb[prefix] = peak_mem
                baseline_memory_kb[prefix] = (info or {}).get('baseline_memory', 0)

    if not all_data:
        print("  No data was parsed. No CSV will be generated.")
        return

    # --- Merge the data ---
    try:
    
        COMMON_MERGE_COLS = ['Step', 'Haplotype_Name', 'Nodes', 'Edges']
        
        # 1. Define time metric structures for use in merging and filtering (FIX: SCOPE)
        BASE_TIME_METRICS = ['UpdateDS_Time', 'DFS_Time', 'TopologicalOrderTime', 'MinI_Time', 'MaxO_Time', 'SB_Compute_Time', 'Total_Time']
        BASE_VERTEX_METRICS = ['AffectedVertices', 'ExploredVertices', 'TotalAffectedEdges', 'TotalExploredEdges']
        PREFIXES = ['static_', 'incRecompute_', 'incNormalBST_', 'incBST_','incBSTBatch_','incMinIBST_','incMinIBSTH1_',
                    'incRecomputeBatch_', 'incBSTBatchH1_', 'incBSTBatchH2_']
        
        # 2. Perform the full merge
        merge_keys = ['Step']
        if mode == "peredge":
            merge_keys = ['Step', 'Haplotype_Name']

        merged_df = all_data[0].copy().set_index(merge_keys)

        for df in all_data[1:]:
            df_to_join = df.set_index(merge_keys)

            # join only metric columns (not the shared identity columns)
            cols_to_join = [col for col in df_to_join.columns
                            if col not in COMMON_MERGE_COLS and col != 'LocalEdgeNo']
            merged_df = merged_df.join(df_to_join[cols_to_join], how='outer')

        merged_df = merged_df.reset_index()

        
        # 3. Collect all time/compute columns
        metric_cols = []

        for prefix in PREFIXES:
            for metric in BASE_TIME_METRICS + BASE_VERTEX_METRICS:
                col = f'{prefix}{metric}'
                if col in merged_df.columns:
                    metric_cols.append(col)

        # 4. Separate the 'Base' step (must be preserved)
        non_base_df = merged_df[merged_df['Step'] != 'Base'].copy()
        base_df = merged_df[merged_df['Step'] == 'Base'].copy()
        
        # 5. Filter non-base rows: a genome/edge without any time is dropped.
        # Only the TIME columns decide (not the vertex/edge counts: a count
        # that overflowed in a program, e.g. a 32-bit sum over a whole HPRC
        # genome, is negative and must not cancel the times).
        valid_time_cols = [col for col in metric_cols
                           if col in non_base_df.columns and col.endswith('_Time')]

        for col in [c for c in metric_cols if c in non_base_df.columns]:
             non_base_df[col] = pd.to_numeric(non_base_df[col], errors='coerce').fillna(0.0)

        # Calculate the sum and filter
        time_sum = non_base_df[valid_time_cols].sum(axis=1)
        filtered_non_base_df = non_base_df[time_sum > 1e-9]

        # 6. Recombine the Base row with the filtered non-base rows
        merged_df = pd.concat([base_df, filtered_non_base_df]).reset_index(drop=True)
        
        print(f"  ℹ️ Filtered out {len(non_base_df) - len(filtered_non_base_df)} haplotypes with zero time.")
        
        
        
        
        CASE_METRICS = ['Case1', 'Case2', 'Case3', 'Case4', 'Case5']
        
        rename_dict = {}
        # the edge-by-edge algorithms classify every edge identically: the
        # cases are taken from the first of them with case columns (incBST;
        # in perbatch incMinIBSTH1), the others' case columns are dropped
        case_source = next((a for a in ('incBST', 'incMinIBST', 'incMinIBSTH1')
                            if f'{a}_Case1' in merged_df.columns), None)
        for case in CASE_METRICS:
            for other in ('incRecompute', 'incNormalBST', 'incBST', 'incMinIBST', 'incMinIBSTH1'):
                if other != case_source and f'{other}_{case}' in merged_df.columns:
                    merged_df = merged_df.drop(columns=[f'{other}_{case}'], errors='ignore')
            # Rename the source's columns to the final CaseX name
            if case_source and f'{case_source}_{case}' in merged_df.columns:
                rename_dict[f'{case_source}_{case}'] = case
        
        merged_df = merged_df.rename(columns=rename_dict)
        
        # --- PascalCase Renaming and Final Column Ordering (CRITICAL FOR NEW METRICS) ---
        
        RENAME_MAP = {
            'Haplotype_Name': 'HaplotypeName',
            'static_SB_Count': 'SuperbubbleCount',
        }

        # Function to convert snake_case column names (e.g., static_DFS_Time) to PascalCase (staticDFSTime)
        def to_pascal_case(name):
            if name in RENAME_MAP:
                return RENAME_MAP[name]

            # Preserve the prefixed parts exactly (so acronyms like DFS, MinI, MaxO stay uppercase)
            # and convert by simply removing the underscores between parts.
            # Special-case the Time_s sentinel which we want as e.g. staticTimeS
            if name.endswith('Time_s'):
                # e.g. 'static_Time_s' -> 'staticTimeS'
                prefix = name.split('_')[0]
                return prefix + 'TimeS'

            if name.startswith(('static_', 'incRecompute_', 'incNormalBST_', 'incBST_', 'incBSTBatch_', 'incMinIBST_', 'incMinIBSTH1_',
                                'incRecomputeBatch_', 'incBSTBatchH1_', 'incBSTBatchH2_')):
                parts = name.split('_')
                # Keep inner parts exactly (no .capitalize()), just glue them:
                # 'static_DFS_Time' -> ['static','DFS','Time'] -> 'static' + 'DFS' + 'Time' -> 'staticDFSTime'
                return parts[0] + ''.join(parts[1:])

            return name  # Case1-5 and common columns

        new_columns = {col: to_pascal_case(col) for col in merged_df.columns}
        merged_df = merged_df.rename(columns=new_columns)
        
        # --- Define Final Order using the new PascalCase names ---
        
        BASE_COLS = ['Step', 'HaplotypeName', 'Nodes', 'Edges']
        TIME_METRICS = ['UpdateDSTime', 'DFSTime', 'TopologicalOrderTime', 'MinITime', 'MaxOTime', 'SBComputeTime', 'TotalTime']
        VERTEX_METRICS = ['AffectedVertices', 'ExploredVertices', 'TotalAffectedEdges', 'TotalExploredEdges']
        GRAPH_METRICS = ['TotalVertices', 'TotalEdges']
        # the batch algorithm writes what a batch did instead of the per-edge cases
        BATCH_METRICS = ['NewVertices', 'SourcesPlaced', 'BackwardEdges', 'Rounds', 'Adoptions']
        PREFIXES = ['static', 'incRecompute', 'incNormalBST', 'incBST', 'incBSTBatch','incMinIBST','incMinIBSTH1',
                    'incRecomputeBatch', 'incBSTBatchH1', 'incBSTBatchH2']

        final_order = BASE_COLS + ['SuperbubbleCount']

        for metric in VERTEX_METRICS:
            for prefix in PREFIXES:
                col_name = f'{prefix}{metric}'
                if col_name in merged_df.columns:
                    final_order.append(col_name)
        
        # 1. Add all detailed time metrics (3 columns each)
        for metric in TIME_METRICS:
             for prefix in PREFIXES:
                col_name = f'{prefix}{metric}'
                if col_name in merged_df.columns:
                    final_order.append(col_name)

        for metric in GRAPH_METRICS:
            for prefix in PREFIXES:
                col_name = f'{prefix}{metric}'
                if col_name in merged_df.columns:
                    final_order.append(col_name)            

        for metric in BATCH_METRICS:
            for prefix in ('incRecomputeBatch', 'incBSTBatch', 'incBSTBatchH1', 'incBSTBatchH2'):
                col_name = f'{prefix}{metric}'
                if col_name in merged_df.columns:
                    final_order.append(col_name)
        
        
        # 3. Add Case metrics (now just the 5 CaseX columns)
        final_order.extend(CASE_METRICS)
        
        final_columns = [col for col in final_order if col in merged_df.columns]

        # Ensure correct column subsetting
        merged_df = merged_df[final_columns]
        
        # Sort logic
        merged_df['SortStep'] = merged_df['Step'].apply(lambda x: 0 if x == 'Base' else (int(x) if pd.notna(x) else float('inf')))
        merged_df = merged_df.sort_values(by='SortStep').drop(columns=['SortStep']).reset_index(drop=True)

        # --- Type and NaN Handling ---
        for col in merged_df.columns:
            if 'Time' in col or 'Compute' in col:
                merged_df[col] = pd.to_numeric(merged_df[col], errors='coerce').astype(float).fillna(0.0)
            elif col in ['SuperbubbleCount', 'Nodes', 'Edges'] or col in CASE_METRICS or any(col.endswith(s) for s in ['AffectedVertices', 'ExploredVertices', 'TotalAffectedEdges', 'TotalExploredEdges'] + BATCH_METRICS):
                merged_df[col] = pd.to_numeric(merged_df[col], errors='coerce').astype('Int64').fillna(0)
            
        merged_df['HaplotypeName'] = merged_df['HaplotypeName'].fillna('N/A')

        # --- Save to CSV ---
        output_filename = time_dir_path / "haplotype_build_comparison.csv"
        merged_df.to_csv(output_filename, index=False)

        # Peak memory per algorithm (PeakMemory(KB) line of each time file),
        # next to the CSV; read by csvToLatex.py (table_PeakMemory.tex)
        memory_filename = time_dir_path / "peak_memory.csv"
        with open(memory_filename, "w") as mf:
            mf.write("Algorithm,PeakMemoryKB,BaselineMemoryKB\n")
            for algo_id, kb in peak_memory_kb.items():
                mf.write(f"{algo_id},{kb},{baseline_memory_kb.get(algo_id, 0)}\n")
        
        # Genomes and final graph size, next to the CSV; read by csvToLatex.py
        # (table_size_haplotypes.tex). Counted from the time file, so genomes
        # that add nothing (dropped from the CSV: zero time) are included.
        # Real*: the input graph (edge list) without the artificial source,
        # sink and their edges; Graph*: the program's graph after the last
        # genome, with them (the denominators of the ratio tables).
        if dataset_info:
            real = real_graph_size(time_dir_path) or ("", "")
            with open(time_dir_path / "dataset_info.csv", "w") as inf:
                inf.write("Genomes,RealNodes,RealEdges,GraphNodes,GraphEdges\n")
                inf.write(f"{dataset_info['genomes']},{real[0]},{real[1]},"
                          f"{dataset_info['final_nodes']},{dataset_info['final_edges']}\n")

        print(f"  ✅ Successfully saved data to: {output_filename}")

    except Exception as e:
        print(f"  ❌ Error merging or saving data for {time_dir_path}: {e}", file=sys.stderr)


# ==============================================================================
#  Main Function 
# ==============================================================================

def main():
    """Finds the timing folders and writes their CSVs; see the guide at the top of this file."""
    parser = argparse.ArgumentParser(
        description="Parses haplotype time-tracking files into CSVs.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-d", "--dataset",
        choices=["ecoli", "human", "hprc", "hprc2", "all"],
        required=True,
        help="The dataset to process. 'ecoli', 'human', 'hprc' (Minigraph-Cactus, "
             "experimental), or 'all' (= ecoli + human)."
    )
    parser.add_argument("-c", "--chromosome",
        nargs='+', 
        help="Specific chromosome numbers (e.g., 1 2 21)."
    )
    parser.add_argument(
        "-m", "--mode",
        choices=["perbatch", "peredge"],
        default="perbatch",
        help="Parsing mode based on file format (counts before or after data)."
    )

    args = parser.parse_args()

    try:
        # 1. Define the name of the repository we are looking for
        REPO_NAME = "IncSuperbubble"
        
        # 2. Get the directory containing this script
        try:
            # Standard way (when running as script)
            script_dir = pathlib.Path(__file__).resolve().parent
        except NameError:
            # Fallback (for interactive, etc.)
            script_dir = pathlib.Path.cwd()

        # 3. Find the repo root path object
        REPO_ROOT_PATH = find_repo_root(script_dir, REPO_NAME)
        OUTPUT_DIR = REPO_ROOT_PATH / "output"
        
        print(f"ℹ️ Located repository root: {REPO_ROOT_PATH}")
        print(f"ℹ️ Looking for output in: {OUTPUT_DIR}")

    except FileNotFoundError as e:
        print(f"❌ Critical Error: {e}", file=sys.stderr)
        print("Please ensure this script is located somewhere inside the 'IncSuperbubble' repository.", file=sys.stderr)
        sys.exit(1)

    # --- 1. Process Ecoli Dataset ---
    if args.dataset in ["ecoli", "all"]:
        ecoli_time_base = OUTPUT_DIR / "Ecoli" / "timeFile"
        if ecoli_time_base.exists():
            # Check the main folder and the specific sub-folders
            sub_folders = [ecoli_time_base, ecoli_time_base / "perEdge", ecoli_time_base / "perbatch"]
            
            for folder in sub_folders:
                if folder.exists():
                    process_time_directory(folder,args.mode)
        else:
            print(f"\nℹ️ 'Ecoli/timeFile' directory not found at {ecoli_time_base}, skipping.")
    """
    if args.dataset in ["ecoli", "all"]:
        ecoli_time_dir = OUTPUT_DIR / "Ecoli" / "timeFile"
        process_time_directory(ecoli_time_dir)
    """
    # --- 2. Process Human Dataset ---
    # (HPRC: Minigraph-Cactus chromosomes, same layout as Human)
    groups = ((["Human"] if args.dataset in ["human", "all"] else []) + (["HPRC"] if args.dataset == "hprc" else [])
              + (["HPRC2"] if args.dataset == "hprc2" else []))
    for group in groups:
        human_output_dir = OUTPUT_DIR / group
        if not human_output_dir.is_dir():
            print(f"\nℹ️ '{group}' directory not found at {human_output_dir}, skipping.")
        else:
            print(f"\nScanning {human_output_dir} for chromosome directories...")
            # Find all directories starting with "chromosome_"
            chromosome_dirs = sorted([d for d in human_output_dir.iterdir() if d.is_dir() and d.name.startswith("chromosome_")])
            
            if args.chromosome:
                # Converts ['1', '2'] into ['chromosome_1', 'chromosome_2']
                target_names = [f"chromosome_{c}" for c in args.chromosome]
                chromosome_dirs = [d for d in chromosome_dirs if d.name in target_names]

            for chrom_dir in chromosome_dirs:
                time_base = chrom_dir / "timeFile"
                if not time_base.exists():
                    continue
                
                # Define exactly which sub-folders to check
                # This ensures we check timeFile/, timeFile/perEdge/, and timeFile/perbatch/
                sub_folders = [time_base, time_base / "perEdge", time_base / "perbatch"]
                
                for folder in sub_folders:
                    if folder.exists():
                        process_time_directory(folder,args.mode)
            """
            if not chromosome_dirs:
                print("  No 'chromosome_*' directories found.")
            else:
                print(f"  Found {len(chromosome_dirs)} chromosome directories to process.")
                
                # Loop over each chromosome directory and process its 'timeFile' subdir
                for chrom_dir in chromosome_dirs:
                    human_time_dir = chrom_dir / "timeFile"
                    process_time_directory(human_time_dir)
            """
    print("\n🎉 All processing complete!")


if __name__ == "__main__":
    main()