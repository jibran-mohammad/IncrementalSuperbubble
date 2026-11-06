#!/usr/bin/env python3
"""
compareBubbleFinder.py - how much the DAG conversion (mcToDag.py) changes the
superbubbles of the HPRC graphs: BubbleFinder on the original GFA (bidirected)
against BubbleFinder on the converted DAG (directed), and our own flag
dag_affected (a site is touched if the conversion removed something at or
inside it) against the real answer (does the same superbubble exist in the
original graph?).

Input, one folder per graph under --dir (made by bf_compare_all.sh on the
server, see the end of this docstring):
    orig_sb.txt, dag_sb.txt     BubbleFinder's lists (first line: the count,
                                then one "entrance exit" pair per line; the
                                DAG list in new node ids with a sign, the
                                original in segment names)
    chromosome_<c>.nodemap.tsv  new id -> original name (mcToDag.py)
    orig.log, dag.log           the BubbleFinder logs (segments and links)
and the lifespans with dag_affected of the graph's history
(output/<dataset>/chromosome_<c>/history/chr<c>_lifespans.tsv).

A pair is compared as an unordered pair of original segment names.

Run from the repository root:
    python3 Preprocessing/compareBubbleFinder.py --dir output/bf_compare_all \\
        --outdir output/paperTables/dag
writes <outdir>/table_dag_vs_original.tex, <outdir>/table_dag_conversion.tex (every
count of the mcToDag.py stats) and prints the numbers.

On the server (run once, about 10 minutes):
    scp bf_compare_all.sh <server>: ; ~/bf_compare_all.sh
"""
import argparse
import csv
import glob
import json
import os
import pathlib
import re
import sys

REPO = pathlib.Path(__file__).resolve().parent.parent

# label of the folder under --dir, column label, history folder, lifespans file, group
GRAPHS = [("chr11", "chr11", "HPRC/chromosome_11", "chr11", "Medium"),
          ("chr17", "chr17", "HPRC/chromosome_17", "chr17", "Medium"),
          ("chr18", "chr18", "HPRC/chromosome_18", "chr18", "Medium"),
          ("chr22", "chr22", "HPRC/chromosome_22", "chr22", "Medium"),
          ("chrX", "chrX", "HPRC/chromosome_X", "chrX", "Medium"),
          ("chr22v2", "chr22 v2", "HPRC2/chromosome_22", "chr22", "Large")]


def read_pairs(path):
    out = []
    with open(path) as f:
        n = int(f.readline())
        for line in f:
            a, b = line.split()
            out.append((a, b))
    assert len(out) == n, (path, len(out), n)
    return out


def log_counts(path):
    """(segments, links) from the GFA parser line of a BubbleFinder log, None if there is none"""
    if not os.path.exists(path):
        return None
    m = re.search(r"GFA parser: (\d+) segments, (\d+) links", open(path).read())
    return (int(m.group(1)), int(m.group(2))) if m else None


def analyse(folder, hist, chrom):
    nodemap = glob.glob(os.path.join(folder, "*.nodemap.tsv"))[0]
    new2orig = {}
    with open(nodemap) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            n = int(row["new"])
            if n >= 0:
                new2orig[n] = row["orig"]
    orig = read_pairs(os.path.join(folder, "orig_sb.txt"))
    dag = read_pairs(os.path.join(folder, "dag_sb.txt"))
    O = {frozenset(p) for p in orig}
    preserved, D = {}, set()
    for a, b in dag:
        na, nb = int(a.rstrip("+-")), int(b.rstrip("+-"))
        pair = frozenset((new2orig[na], new2orig[nb]))
        D.add(pair)
        preserved[(na, nb)] = pair in O
    res = {"orig": len(orig), "dag": len(dag), "both": len(O & D), "only_orig": len(O - D), "only_dag": len(D - O)}
    c_orig, c_dag = log_counts(os.path.join(folder, "orig.log")), log_counts(os.path.join(folder, "dag.log"))
    res["segments_removed"] = c_orig[0] - c_dag[0] if c_orig and c_dag else None
    res["links_removed"] = c_orig[1] - c_dag[1] if c_orig and c_dag else None
    # our sites with an inside vertex, with the dag_affected flag (edge-list id = new id + 2)
    life = REPO / "output" / hist / "history" / f"{chrom}_lifespans.tsv"
    sites = touched = changed = missed = missing = 0
    with open(life) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["died_step"] != "NA":
                continue
            key = (int(row["entrance"]) - 2, int(row["exit"]) - 2)
            if key not in preserved:
                missing += 1
                continue
            sites += 1
            flag = row["dag_affected"] == "1"
            touched += flag
            changed += not preserved[key]
            missed += (not preserved[key]) and not flag
    res.update(sites=sites, touched=touched, changed=changed, missed=missed, missing=missing)
    return res


def num(x):
    return "-" if x is None else str(x)


def conversion_table(cols, outdir):
    """what the conversion did: every count of the <chromosome>.stats.json of mcToDag.py, one column per graph"""
    groups = []
    for _, g, _ in cols:
        if g not in groups:
            groups.append(g)
    spec = "|l|" + "||".join("|".join("c" for _, gg, _ in cols if gg == g) for g in groups) + "|"
    heads = []
    for i, g in enumerate(groups):
        k = sum(1 for _, gg, _ in cols if gg == g)
        heads.append(rf"\multicolumn{{{k}}}{{c{'|' if i == len(groups) - 1 else '||'}}}{{\textbf{{{g}}}}}")
    lines = [r"\begin{table}[ht]", r"\centering", r"\begin{adjustbox}{max width=\linewidth}",
             rf"\begin{{tabular}}{{{spec}}}", r"\hline", " & " + " & ".join(heads) + r" \\", r"\hline",
             r"\textbf{Dataset} & " + " & ".join(rf"\textbf{{{l}}}" for l, _, _ in cols) + r" \\", r"\hline\hline"]
    # (block, [(row name, section, key)]); a key with pct prints a percentage
    blocks = [("Input graph", [("Segments", "input", "nodes"), ("Links", "input", "links"), ("Base pairs", "input", "bp"),
                               ("Haplotype walks", "input", "haplotype_walks"), ("Reference paths", "input", "reference_paths")]),
              ("Orientation", [("Segments flipped", "orientation", "nodes_flipped"),
                               (r"Segments flipped (\%)", "orientation", "pct_nodes_flipped"),
                               ("Walks reversed", "orientation", "walks_reversed"),
                               ("Walks without a reference segment", "orientation", "walks_without_reference_node")]),
              ("Links removed", [("Inversion links", "links_removed", "inversion_links"),
                                 (r"Inversion links (\%)", "links_removed", "pct_inversion"),
                                 ("Self-loops", "links_removed", "self_loops"),
                                 ("Duplicates merged after orientation", "links_removed", "merged_duplicates_after_orientation"),
                                 ("Back edges (cycles)", "links_removed", "cycle_breaking_back_edges"),
                                 ("Back edges on the reference", "links_removed", "back_edges_on_reference"),
                                 ("With small components", "links_removed", "dropped_with_small_components"),
                                 ("Total", "links_removed", "total_removed"),
                                 (r"Total (\%)", "links_removed", "pct_total_removed")]),
              ("Segments removed", [("Components before the filter", "nodes_removed", "components_before_filter"),
                                    ("Segments", "nodes_removed", "nodes"), (r"Segments (\%)", "nodes_removed", "pct_nodes"),
                                    ("Base pairs", "nodes_removed", "bp"), (r"Base pairs (\%)", "nodes_removed", "pct_bp")]),
              ("Haplotype walks", [("Steps", "walks", "transitions"), ("Broken steps", "walks", "broken_transitions"),
                                   (r"Broken steps (\%)", "walks", "pct_broken"),
                                   ("Haplotypes with a broken step", "walks", "haplotypes_with_a_broken_transition"),
                                   ("Batches", "walks", "batches"), ("Links on no walk", "walks", "unwalked_edges")]),
              ("Output DAG", [("Vertices", "output", "nodes"), ("Edges", "output", "edges"), ("Sources", "output", "sources"),
                              ("Sinks", "output", "sinks"), ("Acyclic", "output", "acyclic")])]
    for bi, (name, rows) in enumerate(blocks):
        lines += [rf"\textbf{{{name}}}" + " &" * len(cols) + r" \\", r"\hline"]
        for row, sec, key in rows:
            cells = []
            for _, _, r in cols:
                v = r["stats"][sec][key]
                cells.append("yes" if v is True else "no" if v is False else f"{v:.3f}" if key.startswith("pct") else str(v))
            lines.append(row + " & " + " & ".join(cells) + r" \\")
        lines.append(r"\hline" if bi == len(blocks) - 1 else r"\hline\hline")
    lines += [r"\end{tabular}", r"\end{adjustbox}",
              r"\caption{What the conversion to a DAG did to the HPRC graphs (medium: release 1, large: release 2).}",
              r"\label{tab:dag_conversion}", r"\end{table}"]
    out = os.path.join(outdir, "table_dag_conversion.tex")
    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {out}")


def main():
    ap = argparse.ArgumentParser(description="DAG conversion: BubbleFinder on the original and the DAG")
    ap.add_argument("--dir", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    cols = []                                       # (column label, group, result)
    for folder, label, hist, chrom, group in GRAPHS:
        d = os.path.join(args.dir, folder)
        if not os.path.exists(os.path.join(d, "dag_sb.txt")):
            print(f"[!] {folder}: no results yet, left out")
            continue
        r = analyse(d, hist, chrom)
        r["stats"] = json.load(open(glob.glob(os.path.join(d, "*.stats.json"))[0]))
        cols.append((label, group, r))
        print(f"{label}: original {r['orig']}, DAG {r['dag']}, both {r['both']}, only original {r['only_orig']}, "
              f"only DAG {r['only_dag']}; sites {r['sites']}, touched {r['touched']}, really only in the DAG "
              f"{r['changed']}, untouched but only in the DAG {r['missed']} (not found: {r['missing']})")
    if not cols:
        sys.exit("no results")
    conversion_table(cols, args.outdir)

    groups = []
    for _, g, _ in cols:
        if g not in groups:
            groups.append(g)
    spec = "|l|" + "||".join("|".join("c" for _, gg, _ in cols if gg == g) for g in groups) + "|"
    lines = [r"\begin{table}[ht]", r"\centering", r"\begin{adjustbox}{max width=\linewidth}",
             rf"\begin{{tabular}}{{{spec}}}", r"\hline"]
    heads = []
    for i, g in enumerate(groups):
        k = sum(1 for _, gg, _ in cols if gg == g)
        heads.append(rf"\multicolumn{{{k}}}{{c{'|' if i == len(groups) - 1 else '||'}}}{{\textbf{{{g}}}}}")
    lines += [" & " + " & ".join(heads) + r" \\", r"\hline",
              r"\textbf{Dataset} & " + " & ".join(rf"\textbf{{{l}}}" for l, _, _ in cols) + r" \\", r"\hline"]
    rows = [("Segments removed", "segments_removed"), ("Links removed", "links_removed"), None,
            ("OG superbubbles", "orig"), ("DAG superbubbles", "dag"), ("In both", "both"),
            ("Only in OG", "only_orig"), ("Only in DAG", "only_dag"), None,
            ("Non-trivial superbubbles", "sites"), ("In OG, not in DAG", "touched"),
            ("Non-trivial, only in DAG", "changed"), ("Changed, not flagged", "missed")]
    for row in rows:
        if row is None:
            lines.append(r"\hline")
            continue
        name, key = row
        lines.append(f"{name} & " + " & ".join(num(r[key]) for _, _, r in cols) + r" \\")
        if key == "only_dag":                       # the share of the original that survives the conversion
            lines.append(r"Preserved (\%) & " + " & ".join(
                f"{100.0 * r['both'] / r['orig']:.2f}" for _, _, r in cols) + r" \\")
    lines += [r"\hline", r"\end{tabular}", r"\end{adjustbox}",
              r"\caption{Superbubbles of the original HPRC graphs (OG, bidirected) and of the converted DAGs "
              r"(directed), found by BubbleFinder, and our non-trivial superbubbles (an inside vertex): "
              r"flagged as touched by the conversion (\texttt{dag\_affected}) and only in the DAG.}",
              r"\label{tab:dag_vs_original}", r"\end{table}"]
    out = os.path.join(args.outdir, "table_dag_vs_original.tex")
    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
