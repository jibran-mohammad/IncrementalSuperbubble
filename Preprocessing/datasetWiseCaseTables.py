import pandas as pd
import os
import argparse

# ============================================================
# CONFIGURATION
# ============================================================

CASE_COLUMNS = ["Case1", "Case2", "Case3", "Case4", "Case5"]

METRICS = {
    "Insert": "InsertTime",
    "DFS": "DFSTime",
    "Topo": "TopologicalOrderTime",
    "MinI": "MinITime",
    "MaxO": "MaxOTime",
    "SBCompute": "SBComputeTime",
}
VERTEX_METRICS = {
    "Affected": ["INCAffectedVertices", "BSTAffectedVertices"],
    "Explored": ["BSTExploredVertices"],
}

ALGORITHMS = ["Simple", "INC", "BST"]


# ============================================================
# Utility Functions
# ============================================================

def compute_speedup(simple_val, other_val):
    """
    Speedup = Simple / Other
    """
    if simple_val == 0 or other_val == 0:
        return float("nan")
    return simple_val / other_val


def compute_percent(simple_val, other_val):
    """
    Percentage reduction relative to Simple
    """
    if simple_val == 0:
        return float("nan")
    return 100.0 * (simple_val - other_val) / simple_val


def format_value(x, precision=6):

    if pd.isna(x):
        return "-"

    if abs(x) < 1e-12:
        x = 0.0

    if x < 0:
        x = 0.0

    if abs(x) < 1e-3:
        return f"{x:.6e}"

    return f"{x:.6f}"


def latex_escape(s):
    return str(s).replace("_", r"\_")


# ============================================================
# Core Logic
# ============================================================

def build_vertex_table(df, metric_name):

    cols = VERTEX_METRICS[metric_name]

    rows = []

    for col in cols:

        row = {"Metric": col.replace("Vertices","").replace("INC","INC ").replace("BST","BST ")}

        for case in CASE_COLUMNS:

            if case not in df.columns:
                row[case] = "-"
                continue

            case_df = df[df[case] > 0]

            if col in case_df.columns:
                row[case] = int(case_df[col].astype(float).sum())
            else:
                row[case] = "-"

        rows.append(row)

    return pd.DataFrame(rows)

    

def build_metric_table(df, metric_name, mode):
    """
    Build one metric table (case-wise).
    Returns DataFrame ready for LaTeX.
    """

    metric_suffix = METRICS[metric_name]

    # Prepare storage:
    # rows_dict[row_name][case] = value
    rows_dict = {}

    for algo in ALGORITHMS:
        rows_dict[f"{algo}"] = {}

    rows_dict["INC_Improvement"] = {}
    rows_dict["BST_Improvement"] = {}

    for case in CASE_COLUMNS:

        if case not in df.columns:
            continue

        case_df = df[df[case] > 0]

        # Compute sums per algorithm
        simple_total = case_df.get(f"Simple{metric_suffix}", pd.Series(dtype=float)).sum()
        inc_total = case_df.get(f"INC{metric_suffix}", pd.Series(dtype=float)).sum()
        bst_total = case_df.get(f"BST{metric_suffix}", pd.Series(dtype=float)).sum()

        # Store raw values
        rows_dict["Simple"][case] = format_value(simple_total)
        rows_dict["INC"][case] = format_value(inc_total)
        rows_dict["BST"][case] = format_value(bst_total)

        if mode == "time":
            inc_imp = compute_speedup(simple_total, inc_total)
            bst_imp = compute_speedup(simple_total, bst_total)
        else:
            inc_imp = compute_percent(simple_total, inc_total)
            bst_imp = compute_percent(simple_total, bst_total)

        rows_dict["INC_Improvement"][case] = format_value(inc_imp, 2)
        rows_dict["BST_Improvement"][case] = format_value(bst_imp, 2)

    # Convert to DataFrame
    final_rows = []
    for row_name, values in rows_dict.items():
        row = {"Metric": row_name}
        for case in CASE_COLUMNS:
            row[case] = values.get(case, "-")
        final_rows.append(row)

    return pd.DataFrame(final_rows)


def generate_latex_table(df, metric_name, dataset_name, out_path):
    """
    Generate LaTeX table file.
    """

    columns = ["Metric"] + CASE_COLUMNS

    colfmt = "|" + "|".join(["c"] * len(columns)) + "|"

    lines = []
    lines.append(r"\begin{table}[ht]")
    lines.append(r"\centering")
    lines.append(r"\begin{tabular}{" + colfmt + r"}")
    lines.append(r"\hline")

    header = [latex_escape(c) for c in columns]
    lines.append(" & ".join([rf"\textbf{{{h}}}" for h in header]) + r" \\")
    lines.append(r"\hline")

    for _, row in df.iterrows():
        row_vals = [latex_escape(str(row[c])) for c in columns]
        lines.append(" & ".join(row_vals) + r" \\")
        lines.append(r"\hline")

    lines.append(r"\end{tabular}")
    lines.append(
        rf"\caption{{{metric_name} case-wise comparison for {dataset_name}}}"
    )
    lines.append(
        rf"\label{{tab:{dataset_name}_{metric_name.lower()}}}"
    )
    lines.append(r"\end{table}")

    with open(out_path, "w") as f:
        f.write("\n".join(lines))

    print(f"✅ Wrote: {out_path}")


# ============================================================
# MAIN
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Generate dataset-wise case tables."
    )

    parser.add_argument(
        "--csv",
        required=True,
        help="Path to haplotype_build_comparison.csv"
    )

    parser.add_argument(
        "--dataset",
        required=True,
        help="Dataset name (e.g., Ecoli, chr11)"
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
        help="time = speedup, percent = percentage improvement"
    )

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.csv)
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="ignore")
    if "Step" in df.columns:
        df = df[df["Step"] != "Base"]

    for metric in METRICS.keys():

        metric_df = build_metric_table(df, metric, args.mode)

        out_file = os.path.join(
            args.outdir,
            f"{args.dataset}_{metric}.tex"
        )

        generate_latex_table(
            metric_df,
            metric,
            args.dataset,
            out_file
        )

    for vmetric in VERTEX_METRICS.keys():

        v_df = build_vertex_table(df, vmetric)

        out_file = os.path.join(
            args.outdir,
            f"{args.dataset}_{vmetric}.tex"
        )

        generate_latex_table(
            v_df,
            vmetric,
            args.dataset,
            out_file
        )    

"""
python3 datasetWiseCaseTables.py \
  --csv /path/to/timeFile/perEdge/haplotype_build_comparison.csv \
  --dataset Ecoli \
  --outdir latex_tables \
  --mode time

or 

python3 datasetWiseCaseTables.py \
  --csv /path/to/timeFile/perbatch/haplotype_build_comparison.csv \
  --dataset Ecoli \
  --outdir latex_tables \
  --mode percent

"""
if __name__ == "__main__":
    main()
