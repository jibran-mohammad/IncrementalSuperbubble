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
    table_<case>_<step>_<mode>.tex   the same per case (Case1..Case5)
    table_cases.tex              number of updates of each case
    table_vertex_ratio.tex, table_edge_ratio.tex
                                 affected / explored vertices and edges
                                 relative to the graph size
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
    Dataset | Case1..Case5 | TotalCases
    """
    rows = []

    for path in csv_paths:
        df = pd.read_csv(path)
        df = normalize_dataframe(df, path)
        name = infer_dataset_name_from_path(path)

        present_cases = [c for c in CASE_COLUMNS_CANDIDATES if c in df.columns]
        if not present_cases:
            # If your cases are stored as algorithm-prefixed, you can extend here later.
            rows.append({"Dataset": name, "TotalCases": "-"})
            continue

        totals = {}
        total_cases_sum = 0
        for c in present_cases:
            s = int(df[c].fillna(0).sum())
            totals[c] = s
            total_cases_sum += s

        row = {"Dataset": name}
        for c in CASE_COLUMNS_CANDIDATES:
            row[c] = str(totals.get(c, "-"))
        row["TotalCases"] = str(total_cases_sum)

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

def build_case_metric_table(csv_paths: list, case_name: str, metric_name: str, mode: str):
    """
    Build table per case per metric.
    """

    rows = []
    metric_suffix = METRICS[metric_name]

    for path in csv_paths:
        df = pd.read_csv(path)
        df = normalize_dataframe(df, path)
        name = infer_dataset_name_from_path(path)

        if case_name not in df.columns:
            continue

        # Filter rows where that case occurred
        case_df = df[df[case_name] > 0]

        if case_df.empty:
            continue

        simple_total = float(case_df[f"static{metric_suffix}"].fillna(0).sum())
        inc_total = float(case_df[f"incRecompute{metric_suffix}"].fillna(0).sum())
        bst_total = float(case_df[f"incBST{metric_suffix}"].fillna(0).sum())
        diffmaxo_total   = float(case_df[f"incMinIBST{metric_suffix}"].fillna(0).sum()) if f"incMinIBST{metric_suffix}" in case_df.columns else float("nan")
        h1_total = float(case_df[f"incMinIBSTH1{metric_suffix}"].fillna(0).sum()) if f"incMinIBSTH1{metric_suffix}" in case_df.columns else float("nan")
        batch_total = float(case_df[f"incBSTBatch{metric_suffix}"].fillna(0).sum()) if f"incBSTBatch{metric_suffix}" in case_df.columns else float("nan")

        inc_speedup = compute_speedup(simple_total, inc_total)
        bst_speedup = compute_speedup(simple_total, bst_total)
        diffmaxo_speedup = compute_speedup(simple_total, diffmaxo_total)

        inc_percent = pct_improvement(simple_total, inc_total)
        bst_percent = pct_improvement(simple_total, bst_total)
        diffmaxo_percent = pct_improvement(simple_total, diffmaxo_total)

        if mode == "time":
            row = {
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
            }
        else:
            row = {
                "Dataset": name,
                "incRecompute_Imp(%)": format_float(inc_percent, precision=2),
                "incBST_Imp(%)": format_float(bst_percent, precision=2),
                "incMinIBST_Imp(%)": format_float(diffmaxo_percent, precision=2),
                "incMinIBSTH1_Imp(%)": format_float(pct_improvement(simple_total, h1_total), precision=2),
                "incBSTBatch_Imp(%)": format_float(pct_improvement(simple_total, batch_total), precision=2),
            }

        rows.append(row)

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

def build_vertex_ratio_table(csv_paths: list):
    rows = []

    for path in csv_paths:
        df = pd.read_csv(path)
        
        # Filter Base BEFORE normalize_dataframe, because perBatch groupby drops Step
        if "Step" in df.columns:
            df = df[df["Step"] != "Base"].copy()
        
        df = normalize_dataframe(df, path)
        name = infer_dataset_name_from_path(path)
        mode_label = detect_mode_from_path(path)

        # No longer need this line since we already filtered above
        # data = df[df["Step"] != "Base"]

        total_verts = int(df["incRecomputeTotalVertices"].fillna(0).sum()) if "incRecomputeTotalVertices" in df.columns else 0
        affected    = int(df["incRecomputeAffectedVertices"].fillna(0).sum()) if "incRecomputeAffectedVertices" in df.columns else 0
        explored    = int(df["incBSTExploredVertices"].fillna(0).sum()) if "incBSTExploredVertices" in df.columns else 0

        affected_ratio = (affected / total_verts * 100) if total_verts > 0 else float("nan")
        explored_ratio = (explored / total_verts * 100) if total_verts > 0 else float("nan")

        rows.append({
            "Dataset":           name,
            "Mode":              mode_label,
            "TotalVertices":     str(total_verts),
            "AffectedVertices":  str(affected),
            "ExploredVertices":  str(explored),
            "Affected/Total(%)": format_float(affected_ratio, precision=2),
            "Explored/Total(%)": format_float(explored_ratio, precision=2),
        })

    return pd.DataFrame(rows) 


def build_edge_ratio_table(csv_paths: list):
    rows = []

    for path in csv_paths:
        df = pd.read_csv(path)
        if "Step" in df.columns:
            df = df[df["Step"] != "Base"].copy()
        df = normalize_dataframe(df, path)
        name = infer_dataset_name_from_path(path)
        mode_label = detect_mode_from_path(path)

        total_edges    = int(df["incRecomputeTotalEdges"].fillna(0).sum())         if "incRecomputeTotalEdges"             in df.columns else 0
        affected_edges = int(df["incRecomputeTotalAffectedEdges"].fillna(0).sum()) if "incRecomputeTotalAffectedEdges"     in df.columns else 0
        explored_edges = int(df["incBSTTotalExploredEdges"].fillna(0).sum()) if "incBSTTotalExploredEdges"     in df.columns else 0

        affected_ratio = (affected_edges / total_edges * 100) if total_edges > 0 else float("nan")
        explored_ratio = (explored_edges / total_edges * 100) if total_edges > 0 else float("nan")

        rows.append({
            "Dataset":              name,
            "Mode":                 mode_label,
            "TotalEdges":           str(total_edges),
            "AffectedEdges":        str(affected_edges),
            "ExploredEdges":        str(explored_edges),
            "Affected/Total(%)":    format_float(affected_ratio, precision=2),
            "Explored/Total(%)":    format_float(explored_ratio, precision=2),
        })

    return pd.DataFrame(rows)        

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

    # -------- cases table --------
    cases_df = build_case_table(csv_paths)
    out_tex = os.path.join(args.outdir, "table_cases.tex")
    generate_latex_table(
        cases_df,
        columns=["Dataset"] + CASE_COLUMNS_CANDIDATES + ["TotalCases"],
        caption="Case distribution across datasets (sum across all steps)",
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

    # -------- case specific metric tables --------
    for case_name in CASE_COLUMNS_CANDIDATES:
        for metric_name in ["Insert","DFS", "TopoOrder", "MinI", "MaxO", "SBCompute"]:

            case_df = build_case_metric_table(csv_paths, case_name, metric_name, mode)

            if case_df.empty:
                continue

            if mode == "time":
                columns = ["Dataset", "static", "incRecompute", "incBST", "incMinIBST", "incMinIBSTH1", "incBSTBatch", "incRecompute_Speedup", "incBST_Speedup", "incMinIBST_Speedup", "incMinIBSTH1_Speedup", "incBSTBatch_Speedup"]
                caption = f"{metric_name} time comparison for {case_name}"
                label = f"tab:{case_name.lower()}_{metric_name.lower()}_time"
            else:
                columns = ["Dataset", "incRecompute_Imp(%)", "incBST_Imp(%)", "incMinIBST_Imp(%)", "incMinIBSTH1_Imp(%)", "incBSTBatch_Imp(%)"]
                caption = f"{metric_name} percentage improvement for {case_name}"
                label = f"tab:{case_name.lower()}_{metric_name.lower()}_percent"

            out_tex = os.path.join(
                args.outdir,
                f"table_{case_name}_{metric_name}_{mode}.tex"
            )

            generate_latex_table(
                case_df,
                columns=columns,
                caption=caption,
                label=label,
                out_tex=out_tex
            )


    # --- Vertex ratio table ---
    vratio_df = build_vertex_ratio_table(csv_paths)
    out_tex = os.path.join(args.outdir, "table_vertex_ratio.tex")
    generate_latex_table(
        vratio_df,
        columns=["Dataset", "Mode", "TotalVertices",
                 "AffectedVertices", "ExploredVertices",
                 "Affected/Total(%)", "Explored/Total(%)"],
        caption="Affected and explored vertices as a percentage of total graph vertices",
        label="tab:vertex_ratio",
        out_tex=out_tex
    )

    # --- Edge ratio table (new) ---
    eratio_df = build_edge_ratio_table(csv_paths)
    out_tex = os.path.join(args.outdir, "table_edge_ratio.tex")
    generate_latex_table(
        eratio_df,
        columns=["Dataset", "Mode", "TotalEdges",
                 "AffectedEdges", "ExploredEdges",
                 "Affected/Total(%)", "Explored/Total(%)"],
        caption="Affected and explored edges as a percentage of total graph edges",
        label="tab:edge_ratio",
        out_tex=out_tex
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
