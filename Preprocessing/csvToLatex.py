r"""
csvToLatex.py

Writes LaTeX summary tables that compare the algorithms over one or more
datasets: total time of each step, speedup over the static algorithm (or
percentage improvement), time per case, and ratios of affected/explored
vertices and edges. One row per dataset.

=====================================================================
HOW TO RUN
=====================================================================

Run from the repository root, with the virtual environment active
(source venv/bin/activate), after Preprocessing/parse_to_csv.py.

    # E. coli and chromosomes 11 and 13, per-edge results, raw times
    python3 Preprocessing/csvToLatex.py \
        --csv output/Ecoli/timeFile/perEdge/haplotype_build_comparison.csv \
              output/Human/chromosome_11/timeFile/perEdge/haplotype_build_comparison.csv \
              output/Human/chromosome_13/timeFile/perEdge/haplotype_build_comparison.csv \
        --outdir output/latexTablesperEdge

    # percentage improvement instead of times and speedups
    python3 Preprocessing/csvToLatex.py --csv <csv files> --outdir <folder> --mode percent

Options:

    --csv      one or more haplotype_build_comparison.csv files; the dataset
               name of each row comes from the path (chromosome_11 -> chr11,
               Ecoli)                                                (required)
    --outdir   folder for the .tex files (created if missing)        (required)
    --mode     time: times and speedup over static;
               percent: improvement over static in %                 (default: time)

Use CSVs of the same mode (all perEdge or all perbatch) in one call.

=====================================================================
OUTPUT (.tex files in --outdir)
=====================================================================

    table_<step>.tex             step = Insert, DFS, TopoOrder, MinI, MaxO,
                                 SBCompute: total time of that step per
                                 algorithm (and speedup / improvement)
    table_TotalTime.tex          total time and peak memory
    table_cases.tex              perEdge: number of updates of each case (and share)
    table_case_time.tex          perEdge: time per case, every algorithm
    table_case_time_best.tex     perEdge: the same, static and the best one
    table_case_steps_<ID>.tex    perEdge: time per case and step, one algorithm
    table_batch.tex              perBatch: what the batch algorithm does
                                 per genome
    table_activity.tex           affected / explored vertices and edges,
                                 total and share of the graph (perEdge:
                                 edge-by-edge; perBatch: edge-by-edge vs
                                 batch)
    table_size_haplotypes.tex    graph size per dataset

Columns use the display names of the algorithms (DISPLAY_NAMES below,
e.g. "Inc-MinI-BST Speedup"); the compared algorithms are listed in
ALGORITHMS below. The comparison with external tools such as BubbleGun is
made by Preprocessing/toolComparisonToLatex.py.
"""

import pandas as pd
import os
import sys
import re

# ----------------------------- CONFIG -----------------------------

# Add algorithms here in future (extensible)
ALGORITHMS = {
    "static": "static",
    "incRecompute": "incRecompute",
    "incBST": "incBST",
    "incMinIBST": "incMinIBST",
    "incMinIBSTH1": "incMinIBSTH1",
    "incBSTBatch": "incBSTBatch"
}

# Names shown in table headers (IDs and meanings: see ALGORITHMS.md)
DISPLAY_NAMES = {
    "static": "Static",
    "incRecompute": "Inc-Recompute",
    "incBST": "Inc-BST",
    "incMinIBST": "Inc-MinI-BST",
    "incMinIBSTH1": "Inc-MinI-BST-H1",
    "incBSTBatch": "Inc-BST-Batch",
}


def display_header(column: str) -> str:
    """'incBST_Speedup' -> 'Inc-BST Speedup'; other columns unchanged."""
    algo, _, rest = column.partition("_")
    if algo in DISPLAY_NAMES:
        return DISPLAY_NAMES[algo] + (" " + rest if rest else "")
    return column

# Metrics -> column suffix used in CSV
METRICS = {
    "Insert": "InsertTime",
    "DFS": "DFSTime",
    "TopoOrder": "TopologicalOrderTime",
    "MinI": "MinITime",
    "MaxO": "MaxOTime",
    "SBCompute": "SBComputeTime",
}

CASE_COLUMNS_CANDIDATES = ["Case1", "Case2", "Case3", "Case4", "Case5"]

PEAK_MEMORY_CANDIDATES = [
    "PeakMemoryMB", "PeakMemory", "Peak_Memory_MB", "Peak_Memory"
]

# ------------------------------------------------------------------


def latex_escape(s: str) -> str:
    """Text for a table cell or header: escapes % & # _ (e.g. "Imp(%)" in the
    percent tables, where an unescaped % would start a comment and break the
    table). $ and braces are left alone: values like $\\approx$123 are LaTeX."""
    s = str(s)
    for ch in ("%", "&", "#", "_"):
        s = s.replace(ch, "\\" + ch)
    return s


def infer_dataset_name_from_path(csv_path: str) -> str:
    parts = csv_path.replace("\\", "/").split("/")

    # 1) Handle chromosome_13, chromosome_1, ...
    for p in parts:
        pl = p.lower()
        m = re.match(r"chromosome_([0-9]+|x|y)$", pl)
        if m:
            return f"chr{m.group(1).upper()}"

    # 2) Handle chr1, chr2, chrX, ...
    for p in parts:
        pl = p.lower()
        if pl.startswith("chr"):
            return p

    # 3) Handle Ecoli
    for p in parts:
        if p.lower() in ["ecoli", "e_coli"]:
            return "Ecoli"

    # fallback: parent folder
    return os.path.basename(os.path.dirname(csv_path))

def find_peak_memory_column(df: pd.DataFrame):
    for col in PEAK_MEMORY_CANDIDATES:
        if col in df.columns:
            return col
    return None


def sum_metric(df: pd.DataFrame, algo_key: str, metric_suffix: str) -> float:
    col = f"{algo_key}{metric_suffix}"
    if col not in df.columns:
        return float("nan")
    return float(df[col].fillna(0).sum())

def compute_speedup(simple_val: float, other_val: float) -> float:
    """
    Speedup = static / Other
    """
    if (
        simple_val == 0
        or other_val == 0
        or pd.isna(simple_val)
        or pd.isna(other_val)
    ):
        return float("nan")
    return simple_val / other_val


def pct_improvement(simple_val: float, other_val: float) -> float:
    # Improvement: how much faster compared to static
    # (static - Other)/static * 100
    if simple_val == 0 or pd.isna(simple_val) or pd.isna(other_val):
        return float("nan")
    return 100.0 * (simple_val - other_val) / simple_val


def format_float(x: float, precision: int = 4) -> str:
    if pd.isna(x):
        return "-"
    return f"{x:.{precision}f}".rstrip('0').rstrip('.')


def generate_latex_table(df: pd.DataFrame, columns: list, caption: str, label: str, out_tex: str):
    """
    Generates normal LaTeX table (not longtable) because these are summary tables.
    """
    # leave out columns that are empty for every dataset (an algorithm that
    # did not run in this mode, e.g. incBSTBatch in perEdge)
    def empty(c):
        return c not in df.columns or all(str(v).strip() in ("-", "nan", "NaN", "") for v in df[c])
    columns = [c for c in columns if c == "Dataset" or not empty(c)]
    num_cols = len(columns)
    colfmt = "|" + "|".join(["c"] * num_cols) + "|"

    lines = []
    lines.append(r"\begin{table}[ht]")
    lines.append(r"\centering")
    lines.append(r"\begin{tabular}{" + colfmt + r"}")
    lines.append(r"\hline")

    header = [latex_escape(display_header(c)) for c in columns]
    lines.append(" & ".join([rf"\textbf{{{h}}}" for h in header]) + r" \\")
    lines.append(r"\hline")

    for _, row in df.iterrows():
        row_vals = []
        for c in columns:
            row_vals.append(latex_escape(row.get(c, "-")))
        lines.append(" & ".join(row_vals) + r" \\")
        lines.append(r"\hline")
    lines.append(r"\end{tabular}")
    lines.append(rf"\caption{{{caption}}}")
    lines.append(rf"\label{{{label}}}")
    lines.append(r"\end{table}")

    with open(out_tex, "w") as f:
        f.write("\n".join(lines))

    print(f"✅ Wrote: {out_tex}")


def build_summary_rows_for_metric(csv_paths: list, metric_name: str, mode: str):
    rows = []
    metric_suffix = METRICS[metric_name]

    for path in csv_paths:
        df = pd.read_csv(path)
        df = normalize_dataframe(df, path)
        name = infer_dataset_name_from_path(path)

        simple_total   = sum_metric(df, ALGORITHMS["static"],   metric_suffix)
        inc_total      = sum_metric(df, ALGORITHMS["incRecompute"],      metric_suffix)
        bst_total      = sum_metric(df, ALGORITHMS["incBST"],      metric_suffix)
        diffmaxo_total = sum_metric(df, ALGORITHMS["incMinIBST"], metric_suffix)
        h1_total = sum_metric(df, ALGORITHMS["incMinIBSTH1"], metric_suffix)
        batch_total = sum_metric(df, ALGORITHMS["incBSTBatch"], metric_suffix)

        inc_speedup      = compute_speedup(simple_total, inc_total)
        bst_speedup      = compute_speedup(simple_total, bst_total)
        diffmaxo_speedup = compute_speedup(simple_total, diffmaxo_total)

        inc_percent      = pct_improvement(simple_total, inc_total)
        bst_percent      = pct_improvement(simple_total, bst_total)
        diffmaxo_percent = pct_improvement(simple_total, diffmaxo_total)

        if mode == "time":
            rows.append({
                "Dataset":          name,
                "static":           format_float(simple_total),
                "incRecompute":              format_float(inc_total),
                "incBST":              format_float(bst_total),
                "incMinIBST":         format_float(diffmaxo_total),
                "incMinIBSTH1": format_float(h1_total),
                "incBSTBatch": format_float(batch_total),
                "incRecompute_Speedup":      format_float(inc_speedup,      precision=2),
                "incBST_Speedup":      format_float(bst_speedup,      precision=2),
                "incMinIBST_Speedup": format_float(diffmaxo_speedup, precision=2),
                "incMinIBSTH1_Speedup": format_float(compute_speedup(simple_total, h1_total), precision=2),
                "incBSTBatch_Speedup": format_float(compute_speedup(simple_total, batch_total), precision=2),
            })
        else:
            rows.append({
                "Dataset":          name,
                "incRecompute_Imp(%)":       format_float(inc_percent,      precision=2),
                "incBST_Imp(%)":       format_float(bst_percent,      precision=2),
                "incMinIBST_Imp(%)":  format_float(diffmaxo_percent, precision=2),
                "incMinIBSTH1_Imp(%)": format_float(pct_improvement(simple_total, h1_total), precision=2),
                "incBSTBatch_Imp(%)": format_float(pct_improvement(simple_total, batch_total), precision=2),
            })

    return pd.DataFrame(rows)

def build_case_table(csv_paths: list):
    """
    Dataset | Edges | Case1..Case5: number of edge updates of each case and
    its share ("289 (3.8%)"). The cases are those of the edge-by-edge
    algorithms (all four classify every edge the same way; the CSV takes
    them from incBST).
    """
    rows = []
    for path in csv_paths:
        df = load_updates(path)
        name = infer_dataset_name_from_path(path)
        present = [c for c in CASE_COLUMNS_CANDIDATES if c in df.columns]
        row = {"Dataset": name}
        if not present:
            row["Edges"] = "-"
            rows.append(row)
            continue
        totals = {c: int(pd.to_numeric(df[c], errors="coerce").fillna(0).sum()) for c in present}
        all_edges = sum(totals.values())
        row["Edges"] = f"{all_edges:,}"
        for c in CASE_COLUMNS_CANDIDATES:
            row[c] = f"{totals[c]:,} ({100.0 * totals[c] / all_edges:.1f}%)" if c in totals and all_edges else "-"
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- cases (perEdge)
# Time per case only exists in perEdge: there every row is one edge with one
# case. (In perBatch a row is a genome with edges of every case.)

CASE_TIME_ALGORITHMS = ["static", "incRecompute", "incBST", "incMinIBST", "incMinIBSTH1"]
BEST_EDGE_BY_EDGE = ["incMinIBSTH1", "incMinIBST", "incBST"]
STEP_COLUMNS = [("Insert", "InsertTime"), ("DFS", "DFSTime"), ("TopoOrder", "TopologicalOrderTime"),
                ("MinI", "MinITime"), ("MaxO", "MaxOTime"), ("SBCompute", "SBComputeTime")]


def load_updates(path):
    """The CSV without the base graph row."""
    df = pd.read_csv(path, low_memory=False)
    if "Step" in df.columns:
        df = df[df["Step"].astype(str) != "Base"]
    return df


def case_of_rows(df):
    """perEdge: the case (1-5) of every row, 0 if a row has none."""
    cols = [c for c in CASE_COLUMNS_CANDIDATES if c in df.columns]
    if not cols:
        return pd.Series(0, index=df.index)
    m = df[cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    case = m.values.argmax(axis=1) + 1
    case[m.values.max(axis=1) <= 0] = 0
    return pd.Series(case, index=df.index)


def num(df, col):
    return pd.to_numeric(df[col], errors="coerce").fillna(0)


def build_case_time_table(csv_paths: list, algorithms: list):
    """
    perEdge: Dataset | Algorithm | Total | Case 1..5: total time (s) of each
    algorithm on the edges of each case, with its share of the algorithm's
    total time and its speedup over static on the same edges, e.g.
    "0.147 (11%, 3.7x)"; Total: "1.388 (10.4x)".
    """
    rows = []

    def speedup(static_t, t):
        return f"{static_t / t:.1f}x" if static_t is not None and t > 0 else None

    for path in csv_paths:
        df = load_updates(path)
        name = infer_dataset_name_from_path(path)
        case = case_of_rows(df)
        static_total, static_case = None, {}
        if "staticTotalTime" in df.columns:
            st = num(df, "staticTotalTime")
            static_total = float(st.sum())
            static_case = {k: float(st[case == k].sum()) for k in range(1, 6)}
        for a in algorithms:
            col = f"{a}TotalTime"
            if col not in df.columns:
                continue
            t = num(df, col)
            total = float(t.sum())
            sp = speedup(static_total, total) if a != "static" else None
            row = {"Dataset": name, "Algorithm": DISPLAY_NAMES.get(a, a),
                   "Total": f"{total:.3f}" + (f" ({sp})" if sp else "")}
            for k in range(1, 6):
                v = float(t[case == k].sum())
                if total <= 0:
                    row[f"Case {k}"] = "-"
                    continue
                spk = speedup(static_case.get(k), v) if a != "static" else None
                row[f"Case {k}"] = f"{v:.3f} ({100.0 * v / total:.0f}%" + (f", {spk})" if spk else ")")
            rows.append(row)
    return pd.DataFrame(rows)


def build_case_step_table(csv_paths: list, algo: str):
    """
    perEdge, one algorithm: Dataset | Case | Insert .. SBCompute | Total: time
    (ms) of each step on the edges of each case. The Total column is the
    algorithm's row of table_case_time (in ms).
    """
    rows = []
    for path in csv_paths:
        df = load_updates(path)
        if f"{algo}TotalTime" not in df.columns:
            continue
        name = infer_dataset_name_from_path(path)
        case = case_of_rows(df)
        for k in range(1, 6):
            sel = df[case == k]
            row = {"Dataset": name, "Case": str(k), "Edges": f"{len(sel):,}"}
            total = 0.0
            for label, suffix in STEP_COLUMNS:
                col = f"{algo}{suffix}"
                v = 1000.0 * float(num(sel, col).sum()) if col in sel.columns else 0.0
                total += v
                row[label] = f"{v:,.1f}"
            row["Total"] = f"{total:,.1f}"
            rows.append(row)
    return pd.DataFrame(rows)


BATCH_TABLE_COLUMNS = [("NewVertices", "New vertices"),
                       ("SourcesPlaced", "Sources placed under their first parent"),
                       ("BackwardEdges", "Backward edges repaired"),
                       ("Rounds", "Reorder rounds"),
                       ("Adoptions", "Adoptions")]


def build_batch_table(csv_paths: list):
    """perBatch: what the batch algorithm did, total over the genomes (mean per genome)."""
    rows = []
    for path in csv_paths:
        df = load_updates(path)
        if "incBSTBatchTotalTime" not in df.columns:
            continue
        name = infer_dataset_name_from_path(path)
        g = len(df)
        row = {"Dataset": name, "Genomes": str(g)}
        for col, label in BATCH_TABLE_COLUMNS:
            c = f"incBSTBatch{col}"
            if c in df.columns and g:
                v = float(num(df, c).sum())
                row[label] = f"{int(v):,} ({v / g:,.1f})"
            else:
                row[label] = "-"
        rows.append(row)
    return pd.DataFrame(rows)

def build_vertex_activity_table(csv_paths):

    """
    Dataset | AffectedVertices | ExploredVertices | ReductionRatio
    """

    rows = []

    for path in csv_paths:
        df = pd.read_csv(path)
        df = normalize_dataframe(df, path)

        name = infer_dataset_name_from_path(path)

        affected = df.get("incRecomputeAffectedVertices", 0).sum()
        explored = df.get("incBSTExploredVertices", 0).sum()

        if explored == 0:
            ratio = float("nan")
        else:
            ratio = affected / explored

        rows.append({
            "Dataset": name,
            "AffectedVertices": format_float(affected, 0),
            "ExploredVertices": format_float(explored, 0),
            "ReductionRatio": format_float(ratio, 2)
        })

    return pd.DataFrame(rows)

def build_total_time_table(csv_paths: list, mode: str):
    """
    Dataset | static_TotalTime | incRecompute_TotalTime | incBST_TotalTime | incRecompute_Imp(%) | incBST_Imp(%) | PeakMemory
    """
    rows = []

    for path in csv_paths:
        df = pd.read_csv(path)
        df = normalize_dataframe(df, path)
        name = infer_dataset_name_from_path(path)

        peak_col = find_peak_memory_column(df)
        peak_mem = float(df[peak_col].max()) if peak_col else float("nan")

        simple_total = sum_metric(df, ALGORITHMS["static"], "TotalTime")
        inc_total = sum_metric(df, ALGORITHMS["incRecompute"], "TotalTime")
        bst_total = sum_metric(df, ALGORITHMS["incBST"], "TotalTime")
        diffmaxo_total   = sum_metric(df, ALGORITHMS["incMinIBST"], "TotalTime")
        h1_total = sum_metric(df, ALGORITHMS["incMinIBSTH1"], "TotalTime")
        batch_total = sum_metric(df, ALGORITHMS["incBSTBatch"], "TotalTime")
        
        diffmaxo_speedup = compute_speedup(simple_total, diffmaxo_total)
        diffmaxo_percent = pct_improvement(simple_total, diffmaxo_total)

        inc_speedup = compute_speedup(simple_total, inc_total)
        bst_speedup = compute_speedup(simple_total, bst_total)

        inc_percent = pct_improvement(simple_total, inc_total)
        bst_percent = pct_improvement(simple_total, bst_total)

        if mode == "time":
            rows.append({
                "Dataset": name,
                "static": format_float(simple_total),
                "incRecompute": format_float(inc_total),
                "incBST": format_float(bst_total),
                "incMinIBST":format_float(diffmaxo_total),
                "incMinIBSTH1": format_float(h1_total),
                "incBSTBatch": format_float(batch_total),
                "incRecompute_Speedup": format_float(inc_speedup, precision=2),
                "incBST_Speedup": format_float(bst_speedup, precision=2),
                "incMinIBST_Speedup": format_float(diffmaxo_speedup, precision=2),
                "incMinIBSTH1_Speedup": format_float(compute_speedup(simple_total, h1_total), precision=2),
                "incBSTBatch_Speedup": format_float(compute_speedup(simple_total, batch_total), precision=2),
                "PeakMemory": format_float(peak_mem, precision=2),
            })
        else:
            rows.append({
                "Dataset": name,
                "incRecompute_Imp(%)": format_float(inc_percent, precision=2),
                "incBST_Imp(%)": format_float(bst_percent, precision=2),
                "incMinIBST_Imp(%)":  format_float(diffmaxo_percent, precision=2),
                "incMinIBSTH1_Imp(%)": format_float(pct_improvement(simple_total, h1_total), precision=2),
                "incBSTBatch_Imp(%)": format_float(pct_improvement(simple_total, batch_total), precision=2),
                "PeakMemory": format_float(peak_mem, precision=2),
            })

    return pd.DataFrame(rows)

def build_memory_table(csv_paths: list):
    """
    Dataset | <algorithm> (MB) ... : peak memory of every algorithm, from
    peak_memory.csv next to each CSV (written by parse_to_csv.py from the
    PeakMemory(KB) line of each <ID>_time.txt)
    """
    rows = []
    for path in csv_paths:
        name = infer_dataset_name_from_path(path)
        mem_path = os.path.join(os.path.dirname(path), "peak_memory.csv")
        row = {"Dataset": name}
        memory = {}
        if os.path.exists(mem_path):
            mem = pd.read_csv(mem_path)
            memory = dict(zip(mem["Algorithm"], mem["PeakMemoryKB"]))
        else:
            print(f"  [!] {name}: no peak_memory.csv (run parse_to_csv.py again)")
        for algo in ALGORITHMS.values():
            kb = memory.get(algo)
            row[algo] = format_float(kb / 1024, precision=1) if kb else "-"
        rows.append(row)
    return pd.DataFrame(rows)


def build_size_haplotype_table(csv_paths: list):
    """
    Dataset | Haplotypes | TotalNodes | TotalEdges: genomes added to the base
    graph and the size of the final graph WITHOUT the artificial source, sink
    and their edges (the input graph), from dataset_info.csv next to each CSV
    (parse_to_csv.py: genomes from the time files, size from the edge list).
    Without that file: from the CSV (genomes that add nothing are missing
    there, and the size includes the artificial vertices and edges).
    """
    rows = []

    for path in csv_paths:
        name = infer_dataset_name_from_path(path)
        info_path = os.path.join(os.path.dirname(path), "dataset_info.csv")
        if os.path.exists(info_path):
            info = pd.read_csv(info_path).iloc[0]
            hap_count = int(info["Genomes"])
            if "RealNodes" in info and pd.notna(info["RealNodes"]):
                total_nodes, total_edges = int(info["RealNodes"]), int(info["RealEdges"])
            else:
                print(f"  [!] {name}: no real graph size in dataset_info.csv; with the artificial source/sink")
                total_nodes, total_edges = int(info["GraphNodes"]), int(info["GraphEdges"])
        else:
            print(f"  [!] {name}: no dataset_info.csv (run parse_to_csv.py again); genomes and size from the CSV")
            df = pd.read_csv(path, low_memory=False)
            genomes = df.loc[df["Step"].astype(str) != "Base", "HaplotypeName"] if "Step" in df.columns else df["HaplotypeName"]
            hap_count = int(genomes.nunique())
            total_nodes = int(df["Nodes"].max()) if "Nodes" in df.columns else 0
            total_edges = int(df["Edges"].max()) if "Edges" in df.columns else 0

        rows.append({
            "Dataset": name,
            "Haplotypes": str(hap_count),
            "TotalNodes": str(total_nodes),
            "TotalEdges": str(total_edges),
        })

    return pd.DataFrame(rows)

def detect_mode_from_path(path):

    p = path.lower()

    if "perbatch" in p:
        return "perBatch"

    if "peredge" in p:
        return "perEdge"

    return "perEdge"

def normalize_dataframe(df, path):

    mode = detect_mode_from_path(path)

    if mode == "perBatch":

        numeric_cols = df.select_dtypes(include="number").columns.difference(["Nodes","Edges"])

        df = (
            df
            .groupby("HaplotypeName")[numeric_cols]
            .sum()
            .reset_index()
        )

    return df   

# the edge-by-edge algorithm shown in the activity table (all of them give the
# same affected / explored counts since they share the case logic)
EDGE_BY_EDGE_FOR_ACTIVITY = ["incMinIBSTH1", "incMinIBST", "incBST"]


def build_activity_table(csv_paths: list):
    """
    How much of the graph the algorithms touch, one row per dataset:
      affected  positions that are reordered (window sizes; for the batch
                algorithm also moved blocks), as vertices and as edges (their
                out-degrees)
      explored  vertices (and their out-degrees) reached by the searches
    each as the mean per update (total / number of updates) and as a share
    of the graph: mean per update / mean graph size (= total / sum of the
    graph size over the updates), e.g. 358 (1.95%) = an update reorders 358
    of the graph's 18,388 vertices (Mean V) on average.
      perEdge   one update per edge; the edge-by-edge algorithm (all four
                give the same counts)
      perBatch  one update per genome; the edge-by-edge algorithm (its counts
                are summed over the genome's edges, so they can exceed 100%:
                the same regions are reordered again and again) next to the
                batch algorithm (once per genome)
    Returns (DataFrame, columns, mode, algorithm IDs).
    """
    mode = detect_mode_from_path(csv_paths[0]) if csv_paths else "perEdge"
    rows = []
    algos_used = []

    def cell(total, graph_sum, updates):
        if total is None or graph_sum <= 0 or updates <= 0:
            return "-"
        mean = total / updates
        mean_text = f"{mean:,.1f}" if mean < 10 else f"{mean:,.0f}"
        return f"{mean_text} ({100.0 * total / graph_sum:.2f}%)"

    for path in csv_paths:
        df = pd.read_csv(path, low_memory=False)
        if "Step" in df.columns:
            df = df[df["Step"].astype(str) != "Base"]
        name = infer_dataset_name_from_path(path)
        algos = [next((a for a in EDGE_BY_EDGE_FOR_ACTIVITY if f"{a}AffectedVertices" in df.columns), None)]
        if mode == "perBatch":
            algos.append("incBSTBatch" if "incBSTBatchAffectedVertices" in df.columns else None)
        size_algo = next((a for a in algos if a), None)
        row = {"Dataset": name, "Updates": str(len(df))}
        if size_algo is None:
            rows.append(row)
            continue
        vsum = float(pd.to_numeric(df[f"{size_algo}TotalVertices"], errors="coerce").fillna(0).sum())
        esum = float(pd.to_numeric(df[f"{size_algo}TotalEdges"], errors="coerce").fillna(0).sum())
        row["Mean V"] = f"{vsum / max(len(df), 1):,.0f}"
        row["Mean E"] = f"{esum / max(len(df), 1):,.0f}"
        updates = len(df)
        for a in algos:
            label = "Batch" if a == "incBSTBatch" else "Edge-by-edge"
            def total(col):
                c = f"{a}{col}" if a else None
                return float(pd.to_numeric(df[c], errors="coerce").fillna(0).sum()) if c in df.columns else None
            row[f"Affected V {label}"] = cell(total("AffectedVertices"), vsum, updates)
            row[f"Affected E {label}"] = cell(total("TotalAffectedEdges"), esum, updates)
            row[f"Explored V {label}"] = cell(total("ExploredVertices"), vsum, updates)
            row[f"Explored E {label}"] = cell(total("TotalExploredEdges"), esum, updates)
            if a and a not in algos_used:
                algos_used.append(a)
        rows.append(row)

    labels = ["Edge-by-edge"] + (["Batch"] if mode == "perBatch" else [])
    for r in rows:
        if mode == "perBatch" and "Updates" in r:
            r["Genomes"] = r.pop("Updates")
    columns = ["Dataset", "Genomes" if mode == "perBatch" else "Updates", "Mean V", "Mean E"]
    for what in ["Affected V", "Affected E", "Explored V", "Explored E"]:
        columns += [f"{what} {l}" for l in labels]
    return pd.DataFrame(rows), columns, mode, algos_used


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--csv",
        nargs="+",
        required=True,
        help="One or more haplotype_build_comparison.csv paths"
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for LaTeX tables"
    )

    parser.add_argument(
        "--mode",
        choices=["time", "percent"],
        default="time",
        help="Show raw times or percentage improvement"
    )

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    csv_paths = args.csv
    mode = args.mode
    # -------- 4 metric tables --------
    for metric_name in ["Insert","DFS", "TopoOrder", "MinI", "MaxO", "SBCompute"]:
        summary_df = build_summary_rows_for_metric(csv_paths, metric_name, mode)

        out_tex = os.path.join(args.outdir, f"table_{metric_name}.tex")
        if mode == "time":
            columns = ["Dataset", "static", "incRecompute", "incBST", "incMinIBST", "incMinIBSTH1", "incBSTBatch", "incRecompute_Speedup", "incBST_Speedup", "incMinIBST_Speedup", "incMinIBSTH1_Speedup", "incBSTBatch_Speedup"] 
            caption = f"{metric_name} total time comparison across datasets"
            label = f"tab:{metric_name.lower()}_summary_time"
        else:
            columns = ["Dataset", "incRecompute_Imp(%)", "incBST_Imp(%)", "incMinIBST_Imp(%)", "incMinIBSTH1_Imp(%)", "incBSTBatch_Imp(%)"] 
            caption = f"{metric_name} percentage improvement across datasets"
            label = f"tab:{metric_name.lower()}_summary_percent"
        generate_latex_table(
            summary_df,
            columns, 
            caption,
            label,
            out_tex=out_tex
        )

    # -------- cases table (perEdge: one case per edge; perBatch has table_batch) --------
    cases_df = build_case_table(csv_paths) if detect_mode_from_path(csv_paths[0]) == "perEdge" else pd.DataFrame()
    out_tex = os.path.join(args.outdir, "table_cases.tex")
    if not cases_df.empty:
        generate_latex_table(
            cases_df,
            columns=["Dataset", "Edges"] + CASE_COLUMNS_CANDIDATES,
            caption=("Edge updates per case (number and share). 1: u new, 2: v new, 3: both new, "
                     "4: u after v (the order is repaired), 5: u before v (the order stays valid). "
                     "All edge-by-edge algorithms classify every edge the same way."),
            label="tab:cases_summary",
            out_tex=out_tex
        )

    # -------- peak memory per algorithm --------
    memory_df = build_memory_table(csv_paths)
    generate_latex_table(
        memory_df,
        columns=["Dataset"] + list(ALGORITHMS.values()),
        caption="Peak memory (MB) of each algorithm (peak resident memory of the whole program, "
                "VmHWM, mean over the runs of superbubbleRunner.py)",
        label="tab:peak_memory",
        out_tex=os.path.join(args.outdir, "table_PeakMemory.tex")
    )

    total_df = build_total_time_table(csv_paths, mode)
    out_tex = os.path.join(args.outdir, "table_TotalTime.tex")
    if mode == "time":
        columns = ["Dataset", "static", "incRecompute", "incBST", "incMinIBST", "incMinIBSTH1", "incBSTBatch", "incRecompute_Speedup", "incBST_Speedup", "incMinIBST_Speedup", "incMinIBSTH1_Speedup", "incBSTBatch_Speedup"]
        caption = "Total time comparison across datasets (sum of the per-update times; peak memory: Table~\\ref{tab:peak_memory})"
        label = "tab:totaltime_summary_time"
    else:
        columns = ["Dataset", "incRecompute_Imp(%)", "incBST_Imp(%)", "incMinIBST_Imp(%)", "incMinIBSTH1_Imp(%)", "incBSTBatch_Imp(%)"]
        caption = "Total time percentage improvement across datasets"
        label = "tab:totaltime_summary_percent"
    generate_latex_table(
        total_df,
        columns,
        caption,
        label,
        out_tex=out_tex
    )

    size_df = build_size_haplotype_table(csv_paths)
    out_tex = os.path.join(args.outdir, "table_size_haplotypes.tex")
    generate_latex_table(
        size_df,
        columns=["Dataset", "Haplotypes", "TotalNodes", "TotalEdges"],
        caption="Number of haplotypes and final graph size (nodes/edges) per dataset",
        label="tab:size_haplotypes",
        out_tex=out_tex
    )

    # -------- cases: time per case and per step (perEdge), batch table (perBatch) --------
    data_mode = detect_mode_from_path(csv_paths[0]) if csv_paths else "perEdge"
    if data_mode == "perEdge":
        generate_latex_table(
            build_case_time_table(csv_paths, CASE_TIME_ALGORITHMS),
            columns=["Dataset", "Algorithm", "Total"] + [f"Case {k}" for k in range(1, 6)],
            caption=("Time per case (s): total time of each algorithm on the edge updates of each "
                     "case; in brackets its share of the algorithm's total time and its speedup "
                     "over static on the same edges (x)."),
            label="tab:case_time",
            out_tex=os.path.join(args.outdir, "table_case_time.tex"))
        generate_latex_table(
            build_case_time_table(csv_paths, ["static"] + BEST_EDGE_BY_EDGE[:1]),
            columns=["Dataset", "Algorithm", "Total"] + [f"Case {k}" for k in range(1, 6)],
            caption=("Time per case (s), static and the best edge-by-edge algorithm: total time on "
                     "the edge updates of each case; in brackets its share of the algorithm's total "
                     "time and its speedup over static on the same edges (x)."),
            label="tab:case_time_best",
            out_tex=os.path.join(args.outdir, "table_case_time_best.tex"))
        for algo in CASE_TIME_ALGORITHMS:
            step_df = build_case_step_table(csv_paths, algo)
            if step_df.empty:
                continue
            generate_latex_table(
                step_df,
                columns=["Dataset", "Case", "Edges"] + [l for l, _ in STEP_COLUMNS] + ["Total"],
                caption=(f"{DISPLAY_NAMES.get(algo, algo)}: time (ms) of each step on the edge updates "
                         "of each case (Total = its row of the time-per-case table)."),
                label=f"tab:case_steps_{algo}",
                out_tex=os.path.join(args.outdir, f"table_case_steps_{algo}.tex"))
    else:
        batch_df = build_batch_table(csv_paths)
        if not batch_df.empty:
            generate_latex_table(
                batch_df,
                columns=["Dataset", "Genomes"] + [l for _, l in BATCH_TABLE_COLUMNS],
                caption=("What the batch algorithm does per genome: total over all genomes and, in "
                         "brackets, the mean per genome. Reorder rounds: one for the backward edges "
                         "of a genome (if it has any) plus one per adoption of an old source."),
                label="tab:batch",
                out_tex=os.path.join(args.outdir, "table_batch.tex"))

    # --- Activity table: affected / explored vertices and edges ---
    act_df, act_cols, act_mode, act_algos = build_activity_table(csv_paths)
    shown = ", ".join(DISPLAY_NAMES.get(x, x) for x in act_algos) or "-"
    if act_mode == "perBatch":
        act_caption = ("Part of the graph touched per genome, on average: affected = reordered "
                       "positions, explored = vertices reached by the searches, as vertices (V) and "
                       "edges (E, their out-degrees); mean per genome and, in brackets, its share of "
                       "the mean graph size (Mean V, Mean E). Edge-by-edge, the work of all edges of a "
                       "genome is added up and can exceed 100\\% (the same regions are reordered "
                       f"repeatedly); the batch algorithm reorders once per genome. Algorithms: {shown}.")
    else:
        act_caption = ("Part of the graph touched per edge update, on average: affected = reordered "
                       "positions, explored = vertices reached by the search, as vertices (V) and "
                       "edges (E, their out-degrees); mean per update and, in brackets, its share of "
                       "the mean graph size (Mean V, Mean E); the edge-by-edge algorithms give the "
                       f"same counts. Algorithm: {shown}.")
    generate_latex_table(
        act_df,
        columns=act_cols,
        caption=act_caption,
        label="tab:activity",
        out_tex=os.path.join(args.outdir, "table_activity.tex")
    )

if __name__ == "__main__":
    main()
    """
    python3 csvToLatex.py 
      --csv 
      /home/dell/Research/Superbubble/Experiment/IncSuperbubble/output/Ecoli/timeFile/haplotype_build_comparison.csv \
      /home/dell/Research/Superbubble/Experiment/IncSuperbubble/output/chr1/timeFile/haplotype_build_comparison.csv \
      /home/dell/Research/Superbubble/Experiment/IncSuperbubble/output/chr2/timeFile/haplotype_build_comparison.csv \
      ... 
      /home/dell/Research/Superbubble/Experiment/IncSuperbubble/output/chr22/timeFile/haplotype_build_comparison.csv \
      --outdir latex_tables  --mode time or --mode percent

    """
