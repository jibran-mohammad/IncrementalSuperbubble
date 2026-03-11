import pandas as pd
import os
import sys
import re

# ----------------------------- CONFIG -----------------------------

# Add algorithms here in future (extensible)
ALGORITHMS = {
    "Simple": "Simple",
    "INC": "INC",
    "BST": "BST"
}

# Metrics -> column suffix used in CSV
METRICS = {
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
    return str(s).replace("_", r"\_")


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
    Speedup = Simple / Other
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
    # Improvement: how much faster compared to Simple
    # (Simple - Other)/Simple * 100
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
    num_cols = len(columns)
    colfmt = "|" + "|".join(["c"] * num_cols) + "|"

    lines = []
    lines.append(r"\begin{table}[ht]")
    lines.append(r"\centering")
    lines.append(r"\begin{tabular}{" + colfmt + r"}")
    lines.append(r"\hline")

    header = [latex_escape(c) for c in columns]
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
    """
    Build rows:
    Dataset | Simple_Total | INC_Total | BST_Total | INC_Imp(%) | BST_Imp(%)
    """
    rows = []
    metric_suffix = METRICS[metric_name]

    for path in csv_paths:
        df = pd.read_csv(path)
        df = normalize_dataframe(df, path)
        name = infer_dataset_name_from_path(path)

        simple_total = sum_metric(df, ALGORITHMS["Simple"], metric_suffix)
        inc_total = sum_metric(df, ALGORITHMS["INC"], metric_suffix)
        bst_total = sum_metric(df, ALGORITHMS["BST"], metric_suffix)

        inc_speedup = compute_speedup(simple_total, inc_total)
        bst_speedup = compute_speedup(simple_total, bst_total)

        inc_percent = pct_improvement(simple_total, inc_total)
        bst_percent = pct_improvement(simple_total, bst_total)
        if mode == "time":
            rows.append({
                "Dataset": name,
                "Simple": format_float(simple_total),
                "INC": format_float(inc_total),
                "BST": format_float(bst_total),
                "INC_Speedup": format_float(inc_speedup, precision=2),
                "BST_Speedup": format_float(bst_speedup, precision=2),
            })
        else:
            rows.append({
                "Dataset": name,
                "INC_Imp(%)": format_float(inc_percent, precision=2),
                "BST_Imp(%)": format_float(bst_percent, precision=2),
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

        affected = df.get("INCAffectedVertices", 0).sum()
        explored = df.get("BSTExploredVertices", 0).sum()

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
    Dataset | Simple_TotalTime | INC_TotalTime | BST_TotalTime | INC_Imp(%) | BST_Imp(%) | PeakMemory
    """
    rows = []

    for path in csv_paths:
        df = pd.read_csv(path)
        df = normalize_dataframe(df, path)
        name = infer_dataset_name_from_path(path)

        peak_col = find_peak_memory_column(df)
        peak_mem = float(df[peak_col].max()) if peak_col else float("nan")

        simple_total = sum_metric(df, ALGORITHMS["Simple"], "TotalTime")
        inc_total = sum_metric(df, ALGORITHMS["INC"], "TotalTime")
        bst_total = sum_metric(df, ALGORITHMS["BST"], "TotalTime")

        inc_speedup = compute_speedup(simple_total, inc_total)
        bst_speedup = compute_speedup(simple_total, bst_total)

        inc_percent = pct_improvement(simple_total, inc_total)
        bst_percent = pct_improvement(simple_total, bst_total)

        if mode == "time":
            rows.append({
                "Dataset": name,
                "Simple": format_float(simple_total),
                "INC": format_float(inc_total),
                "BST": format_float(bst_total),
                "INC_Speedup": format_float(inc_speedup, precision=2),
                "BST_Speedup": format_float(bst_speedup, precision=2),
                "PeakMemory": format_float(peak_mem, precision=2),
            })
        else:
            rows.append({
                "Dataset": name,
                "INC_Imp(%)": format_float(inc_percent, precision=2),
                "BST_Imp(%)": format_float(bst_percent, precision=2),
                "PeakMemory": format_float(peak_mem, precision=2),
            })

    return pd.DataFrame(rows)

def build_size_haplotype_table(csv_paths: list):
    """
    Dataset | Haplotypes | TotalNodes | TotalEdges
    """
    rows = []

    for path in csv_paths:
        df = pd.read_csv(path)
        name = infer_dataset_name_from_path(path)

        # haplotypes
        if "HaplotypeName" in df.columns:
            hap_count = int(df["HaplotypeName"].nunique()) if "HaplotypeName" in df.columns else len(df)
        else:
            hap_count = 0

        # graph size (take max because repeated per step)
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

        simple_total = float(case_df[f"Simple{metric_suffix}"].fillna(0).sum())
        inc_total = float(case_df[f"INC{metric_suffix}"].fillna(0).sum())
        bst_total = float(case_df[f"BST{metric_suffix}"].fillna(0).sum())

        inc_speedup = compute_speedup(simple_total, inc_total)
        bst_speedup = compute_speedup(simple_total, bst_total)

        inc_percent = pct_improvement(simple_total, inc_total)
        bst_percent = pct_improvement(simple_total, bst_total)

        if mode == "time":
            row = {
                "Dataset": name,
                "Simple": format_float(simple_total),
                "INC": format_float(inc_total),
                "BST": format_float(bst_total),
                "INC_Speedup": format_float(inc_speedup, precision=2),
                "BST_Speedup": format_float(bst_speedup, precision=2),
            }
        else:
            row = {
                "Dataset": name,
                "INC_Imp(%)": format_float(inc_percent, precision=2),
                "BST_Imp(%)": format_float(bst_percent, precision=2),
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
    for metric_name in ["DFS", "TopoOrder", "MinI", "MaxO", "SBCompute"]:
        summary_df = build_summary_rows_for_metric(csv_paths, metric_name, mode)

        out_tex = os.path.join(args.outdir, f"table_{metric_name}.tex")
        if mode == "time":
            columns = ["Dataset", "Simple", "INC", "BST","INC_Speedup","BST_Speedup",]
            caption = f"{metric_name} total time comparison across datasets"
            label = f"tab:{metric_name.lower()}_summary_time"
        else:
            columns = ["Dataset", "INC_Imp(%)", "BST_Imp(%)"]
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

    total_df = build_total_time_table(csv_paths, mode)
    out_tex = os.path.join(args.outdir, "table_TotalTime.tex")
    if mode == "time":
        columns = ["Dataset", "Simple", "INC", "BST","INC_Speedup","BST_Speedup", "PeakMemory"]
        caption = "Total time comparison across datasets (with peak memory)"
        label = "tab:totaltime_summary_time"
    else:
        columns = ["Dataset", "INC_Imp(%)", "BST_Imp(%)", "PeakMemory"]
        caption = "Total time percentage improvement across datasets (with peak memory)"
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
        for metric_name in ["DFS", "TopoOrder", "MinI", "MaxO", "SBCompute"]:

            case_df = build_case_metric_table(csv_paths, case_name, metric_name, mode)

            if case_df.empty:
                continue

            if mode == "time":
                columns = ["Dataset", "Simple", "INC", "BST","INC_Speedup","BST_Speedup",]
                caption = f"{metric_name} time comparison for {case_name}"
                label = f"tab:{case_name.lower()}_{metric_name.lower()}_time"
            else:
                columns = ["Dataset", "INC_Imp(%)", "BST_Imp(%)"]
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


    vertex_df = build_vertex_activity_table(csv_paths)

    out_tex = os.path.join(args.outdir, "table_vertex_activity.tex")

    generate_latex_table(
        vertex_df,
        columns=["Dataset","AffectedVertices","ExploredVertices","ReductionRatio"],
        caption="Affected vertices vs explored vertices and their reduction ratio across datasets",
        label="tab:vertex_activity",
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
