#!/usr/bin/env python3
"""
datasetTable.py - the dataset table of the paper: haplotypes and the size of
the final graph (nodes, edges) per dataset, the datasets in three groups
(small: the human chromosomes, medium: HPRC release 1, large: HPRC release 2).
Same layout as the table of the paper (rows Haplotypes / Total Nodes / Total
Edges, one column per dataset), the columns grouped under Small / Medium /
Large; wrapped in an adjustbox so it never gets wider than the text.

Sizes: dataset_info.csv next to the time files of the dataset (parse_to_csv.py:
haplotypes, and the size of the final graph without the artificial source and
sink); without a real size there (no edge list on this machine) the history
summary of the dataset (new_vertices / new_edges summed over the steps).

Run from the repository root:
    python3 Preprocessing/datasetTable.py --outdir output/paperTables/datasets
writes <outdir>/table_datasets.tex
"""
import argparse
import os
import pathlib
import sys

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))
from Preprocessing.paperTables import history_graph_size                     # noqa: E402

CHROMS = ["11", "17", "18", "22", "X"]
GROUPS = [("Small", [(f"chr{c}", f"Human/chromosome_{c}") for c in CHROMS]),
          ("Medium", [(f"chr{c}", f"HPRC/chromosome_{c}") for c in CHROMS]),
          ("Large", [("chr22 v2", "HPRC2/chromosome_22")])]


def size_of(dataset):
    """(haplotypes, nodes, edges) of a dataset, None if unknown"""
    info = REPO / "output" / dataset / "timeFile" / "perbatch" / "dataset_info.csv"
    if not info.exists():
        return None
    row = pd.read_csv(info).iloc[0]
    haps = int(row["Genomes"])
    if pd.notna(row.get("RealNodes")):
        return haps, int(row["RealNodes"]), int(row["RealEdges"])
    real = history_graph_size(str(info))
    return (haps,) + real if real else None


def main():
    ap = argparse.ArgumentParser(description="Dataset table: small / medium / large")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    cols = []                                       # (group, label, (haplotypes, nodes, edges))
    for group, members in GROUPS:
        for label, dataset in members:
            sz = size_of(dataset)
            if sz is None:
                print(f"[!] no size for {dataset}")
            cols.append((group, label, sz))

    n_small = sum(1 for g, _, _ in cols if g == "Small")
    n_med = sum(1 for g, _, _ in cols if g == "Medium")
    n_large = sum(1 for g, _, _ in cols if g == "Large")
    spec = "|l|" + "|".join(["c"] * n_small) + "||" + "|".join(["c"] * n_med) + "||" + "|".join(["c"] * n_large) + "|"

    def cell(sz, i):
        return str(sz[i]) if sz else "-"

    lines = [r"\begin{table}[t]", r"\centering", r"\setlength{\tabcolsep}{8pt}",
             r"\begin{adjustbox}{max width=\linewidth}",
             rf"\begin{{tabular}}{{{spec}}}", r"\hline",
             rf" & \multicolumn{{{n_small}}}{{c||}}{{\textbf{{Small}}}} & \multicolumn{{{n_med}}}{{c||}}"
             rf"{{\textbf{{Medium}}}} & \multicolumn{{{n_large}}}{{c|}}{{\textbf{{Large}}}} \\", r"\hline",
             r"\textbf{Dataset} & " + " & ".join(rf"\textbf{{{l}}}" for _, l, _ in cols) + r" \\", r"\hline"]
    for name, i in (("Haplotypes", 0), ("Total Nodes", 1), ("Total Edges", 2)):
        lines += [rf"\textbf{{{name}}} & " + " & ".join(cell(sz, i) for _, _, sz in cols) + r" \\", r"\hline"]
    lines += [r"\end{tabular}", r"\end{adjustbox}",
              r"\caption{Number of haplotypes and final graph size (nodes/edges, without the artificial source "
              r"and sink) per dataset. Small: human chromosome graphs; medium: the same chromosomes of the HPRC "
              r"release 1 graph; large: chr22 of the HPRC release 2 graph (chr22 v2).}",
              r"\label{tab:datasets}", r"\end{table}"]
    out = os.path.join(args.outdir, "table_datasets.tex")
    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
