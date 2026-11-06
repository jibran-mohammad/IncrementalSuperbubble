#!/usr/bin/env python3
"""
historyTables.py - the paper tables of the superbubble history: what the
genomes do to the superbubbles, the variant types of the final superbubbles,
and (per-edge history) what the edges do, per case of the incremental DFS.

=====================================================================
HOW TO RUN
=====================================================================

Run from the repository root, after superbubbleHistory.py,
superbubbleAlleles.py and (for the per-edge table) edgeHistory.py:

    python3 Preprocessing/historyTables.py \\
        --history output/HPRC/chromosome_11/history/chr11 output/HPRC/chromosome_22/history/chr22 \\
        --label   chr11 chr22 \\
        --edges   output/HPRC/chromosome_22/history/chr22 \\
        --edge-label chr22 \\
        --outdir  output/paperTables/medium/history --tag medium

Options:
    --history     output prefix of superbubbleHistory.py per dataset (reads
                  <prefix>_summary.tsv, <prefix>_lifespans.tsv and, if
                  present, <prefix>_alleles.tsv of superbubbleAlleles.py)
    --label       a column name per dataset
    --edges       output prefix of edgeHistory.py per dataset (reads
                  <prefix>_edge_summary.tsv)
    --edge-label  a column name per --edges dataset
    --edge-note   a sentence added to the caption of the per-edge table,
                  e.g. "Release 2: the first 1.2 million edges."
    --outdir      output folder
    --tag         label suffix (tab:paper_history_<tag>, ...)

=====================================================================
OUTPUT
=====================================================================

    table_history_summary.tex    per dataset: genomes, superbubbles of the
                                 base graph and of the final graph, born,
                                 promoted, died over all genomes, and the
                                 share of the final superbubbles touched by
                                 the DAG normalization (dag_affected)
    table_history_var_types.tex  per dataset: variant type of the final
                                 superbubbles (% of them) and the share of
                                 multi-allelic sites (n_alleles > 2)
    table_edge_history.tex       per dataset and case: edges, edges that
                                 change the superbubbles, superbubbles
                                 destroyed and created, full recomputations
"""
import argparse
import csv
import os
import pathlib
import sys

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))
from Preprocessing.paperTables import write_table                     # noqa: E402

VAR_ORDER = ["SNP", "MNP", "insertion", "deletion", "complex", "mixed", "nested", "nonref", "none"]
CASE_TEXT = {"0": "duplicate", "1": "1", "2": "2", "3": "3", "4": "4", "5": "5", "all": "all"}


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def num(x):
    """thousands separator in LaTeX"""
    return f"{int(x):,}".replace(",", "{,}")


def pct(a, b):
    return f"{100.0 * a / b:.1f}" if b else "-"


def history_summary(prefix):
    """genomes, base, final, born, promoted, died, dag_affected % of the final superbubbles"""
    s = read_tsv(prefix + "_summary.tsv")
    last = s[-1]
    final = [r for r in read_tsv(prefix + "_lifespans.tsv") if r["died_step"] == "NA"]
    dag = [r for r in final if r.get("dag_affected", "NA") not in ("NA", "")]
    affected = sum(1 for r in dag if r["dag_affected"] == "1")
    return {"genomes": len(s) - 1, "base": int(s[0]["alive"]), "final": int(last["alive"]),
            "born": int(last["cum_born"]), "promoted": int(last["cum_promoted"]), "died": int(last["cum_died"]),
            "dag": pct(affected, len(dag)) if dag else "-"}


def var_types(prefix):
    """% of the final superbubbles per variant type, % multi-allelic; None
    without alleles. Only the superbubbles alive at the end of the lifespans
    (those with an inside vertex; the alleles may also list trivial ones)"""
    path = prefix + "_alleles.tsv"
    if not os.path.exists(path):
        return None
    final = {(r["entrance"], r["exit"]) for r in read_tsv(prefix + "_lifespans.tsv") if r["died_step"] == "NA"}
    rows = [r for r in read_tsv(path) if (r["entrance"], r["exit"]) in final]
    out = {t: pct(sum(1 for r in rows if r["var_type"] == t), len(rows)) for t in VAR_ORDER}
    out["multi"] = pct(sum(1 for r in rows if int(r["n_alleles"]) > 2), len(rows))
    out["n"] = len(rows)
    return out


def main():
    ap = argparse.ArgumentParser(description="Paper tables of the superbubble history")
    ap.add_argument("--history", nargs="+", default=[])
    ap.add_argument("--label", nargs="+", default=[])
    ap.add_argument("--edges", nargs="+", default=[])
    ap.add_argument("--edge-label", nargs="+", default=[])
    ap.add_argument("--edge-note", default="")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    tag = f"_{args.tag}" if args.tag else ""
    labels = args.label or [os.path.basename(p) for p in args.history]

    # 1. what the genomes do
    if args.history:
        h = [history_summary(p) for p in args.history]
        rows = [[name] + [cell(x) for x in h] for name, cell in [
            ("Haplotypes", lambda x: num(x["genomes"])),
            ("Base graph", lambda x: num(x["base"])),
            ("Born", lambda x: num(x["born"])),
            ("Promoted", lambda x: num(x["promoted"])),
            ("Died", lambda x: num(x["died"])),
            ("Final graph", lambda x: num(x["final"])),
            ("DAG-affected (\\%)", lambda x: x["dag"])]]
        write_table(os.path.join(args.outdir, "table_history_summary.tex"), ["Dataset"] + labels, rows,
                    "Superbubble history per haplotype (batch algorithm): superbubbles of the base graph, "
                    "born (new sites: containing a vertex the haplotype added), promoted (an existing region "
                    "became a superbubble) and died over all haplotypes, superbubbles of the final graph, and "
                    "the share of the final superbubbles touched by the DAG normalization (an endpoint of a "
                    "removed link inside the site or at its boundary). Base + born + promoted - died = final. "
                    "Only superbubbles with at least one inside vertex are counted; trivial ones (a single "
                    "edge, nothing inside, e.g. every edge of the reference path of the base graph) are not "
                    "variant sites.",
                    f"tab:paper_history{tag}")

        # 2. variant types of the final superbubbles
        v = [var_types(p) for p in args.history]
        if any(x is not None for x in v):
            def vcell(x, k):
                return "-" if x is None else (num(x[k]) if k == "n" else x[k])
            rows = [[t] + [vcell(x, t) for x in v] for t in VAR_ORDER]
            rows += [None, ["Multi-allelic"] + [vcell(x, "multi") for x in v],
                     ["Superbubbles"] + [vcell(x, "n") for x in v]]
            write_table(os.path.join(args.outdir, "table_history_var_types.tex"), ["Variant type (\\%)"] + labels,
                        rows,
                        "Variant type of the final superbubbles (\\% of them), from the haplotype walks: every "
                        "allele (the sequence between entrance and exit) compared with the reference allele "
                        "(SNP, MNP: same length; insertion, deletion: one side empty; complex: different "
                        "lengths); mixed = the non-reference alleles disagree; nested = contains other "
                        "superbubbles; nonref = the reference does not cross the site; none = every haplotype "
                        "crossing the site carries the reference allele. Multi-allelic: more than two distinct "
                        "alleles. Only superbubbles with at least one inside vertex.", f"tab:paper_history_var_types{tag}")

    # 3. what the edges do (per-edge history)
    if args.edges:
        elabels = args.edge_label or [os.path.basename(p) for p in args.edges]
        summ = [{r["case"]: r for r in read_tsv(p + "_edge_summary.tsv")} for p in args.edges]
        cases = [c for c in ["1", "2", "3", "4", "5", "0", "all"] if any(c in s for s in summ)]
        header = ["Case"]
        for lab in elabels:
            header += [f"{lab}: edges", "changing (\\%)", "destroyed", "created", "full"]
        rows = []
        for c in cases:
            row = [CASE_TEXT[c] if c != "all" else "\\textbf{all}"]
            for s in summ:
                r = s.get(c)
                row += (["-"] * 5 if r is None else
                        [num(r["edges"]), pct(int(r["changing"]), int(r["edges"])), num(r["destroyed"]),
                         num(r["created"]), num(r["full"])])
            if c == "all":
                rows.append(None)
            rows.append(row)
        write_table(os.path.join(args.outdir, "table_edge_history.tex"), header, rows,
                    "Superbubble history per inserted edge (H1), by case of the incremental DFS: edges, "
                    "share of the edges that change the superbubbles, superbubbles destroyed and created "
                    "(net change of the edge), and edges for which the superbubbles were recomputed in full "
                    "(u a sink or v a source before the edge). Cases: 1 u new, 2 v new, 3 both new, 4 the "
                    "edge goes against the order (the window between v and u is reordered), 5 the edge "
                    "agrees with the order (no reordering)." + (" " + args.edge_note if args.edge_note else ""),
                    f"tab:paper_edge_history{tag}")


if __name__ == "__main__":
    main()
