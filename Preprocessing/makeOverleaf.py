"""
makeOverleaf.py - builds the Overleaf folder with every table of the paper.

    output/overleaf_tables/tables/
        paper_tables.tex       catalogue of the paper tables: one \\input line
                               each (with label and content), to copy into
                               the text where the table belongs
        appendix_tables.tex    the appendix with ALL tables: first the paper
                               tables (to be moved into the text one by one),
                               then the detailed tables;
                               include with  \\clearpage \\appendix \\input{tables/appendix_tables}
        paper/<mode>/          tables of Preprocessing/paperTables.py
        <mode>/                detailed tables of Preprocessing/csvToLatex.py
        toolComparison/<tool>/ tables of Preprocessing/toolComparisonToLatex.py
    output/overleaf_tables.zip the same folder, to upload

<mode> = <class><Mode>: the class of datasets (small = human chr11, 17, 18,
22, X; medium = HPRC release 1, the same chromosomes; large = HPRC release 2
chr22) and perEdge, perBatch or History (tables of historyTables.py).
Every tabular is wrapped in an adjustbox (never wider than the text; the
large full-time table also not higher than the page), so the preamble needs
\\usepackage{adjustbox}. Labels of the detailed tables get the mode as
prefix (tab:perEdge:..., tab:perBatchHPRC:...) because the same table
exists in several modes; the paper tables keep their own unique labels.

Run from the repository root after the tables are generated:
    python3 Preprocessing/makeOverleaf.py
and upload output/overleaf_tables (or the zip) to Overleaf, replacing the
tables/ folder there.
"""
import os
import pathlib
import re
import shutil

from PIL import Image

REPO = pathlib.Path(__file__).resolve().parent.parent
OUT = REPO / "output"
DEST = OUT / "overleaf_tables"
TABLES = DEST / "tables"

# only the paper tables go into the Overleaf folder. The detailed tables
# (csvToLatex.py) and the old tool comparison tables are left out; to
# include them again, list them here as
#   DETAILED = [("perEdge", OUT / "latexTablesperEdge"),
#               ("perBatch", OUT / "latexTablesperBatch"),
#               ("perBatchHPRC", OUT / "latexTablesperBatch_HPRC")]
#   TOOL_DIRS = ["bubblegun", "bubblefinder", "bubblegun_hprc22", "bubblefinder_hprc22",
#                "bubblegun_hprc22_H1", "bubblefinder_hprc22_H1"]
DETAILED = []
# the three classes of datasets, each with its own tables and plots
CLASSES = [("small", "small datasets: human chr11, chr17, chr18, chr22, chrX"),
           ("medium", "medium datasets: HPRC release 1, chr11, chr17, chr18, chr22, chrX"),
           ("large", "large dataset: HPRC release 2, chr22")]
PAPER = [(f"{c}{m[0].upper()}{m[1:]}", OUT / "paperTables" / c / m)
         for c, _ in CLASSES for m in ("perEdge", "perBatch", "history")]
# the dataset table (haplotypes, nodes, edges; small / medium / large) comes first
PAPER.insert(0, ("datasets", OUT / "paperTables" / "datasets"))
# the long time table in three dataset groups (small / medium / large), perEdge and perBatch;
# it replaces the long tables of the classes (table_full_time)
PAPER.insert(1, ("grouped", OUT / "paperTables" / "grouped"))
# the DAG conversion: BubbleFinder on the original graphs and on the DAGs (HPRC)
PAPER.insert(2, ("dag", OUT / "paperTables" / "dag"))
# medium perEdge: the long table only for chr22, where (nearly) every algorithm finished
PAPER.insert(PAPER.index(("mediumPerEdge", OUT / "paperTables" / "medium" / "perEdge")) + 1,
             ("mediumPerEdgeChr22", OUT / "paperTables" / "medium" / "perEdge_chr22"))
TOOL_DIRS = []
MODE_TITLE = {f"{c}{m[0].upper()}{m[1:]}": f"{text}; {m}" for c, text in CLASSES
              for m in ("perEdge", "perBatch", "history")}
MODE_TITLE["grouped"] = "long time tables: small / medium / large datasets in one table, perEdge and perBatch"
MODE_TITLE["dag"] = "the DAG conversion: superbubbles of the original graphs and of the converted DAGs"
MODE_TITLE["datasets"] = "datasets: haplotypes and graph size, small / medium / large"
MODE_TITLE["mediumPerEdgeChr22"] = "medium datasets: HPRC release 1, chr22 only; perEdge"
# history figures (output/plots/history/<class>/), by the end of the file name
HISTORY_TEXT = {
    "growth": "Superbubbles alive, cumulative born (new sites: containing a vertex the haplotype added), "
              "cumulative promoted (an existing region became a superbubble) and cumulative died, against "
              "the number of haplotypes added; only superbubbles with at least one inside vertex (the base "
              "graph, the reference path, has none)",
    # new_sites (new sites per genome by superpopulation) is left out: the genomes are added in
    # population blocks, so the comparison mostly reflects the position in the order
    "var_types": "Variant type of the final superbubbles, untouched vs touched by the DAG normalization, as "
                 "shares of each group (number of superbubbles above each bar)",
    # edge_destroyed (superbubbles destroyed per inserted edge along the reference) is left out:
    # nearly every edge destroys one, the plot has no visible pattern; the numbers are in
    # table_edge_history
}

# plots in the appendix: only the comparison plots (output/plots/comparison,
# all datasets in one plot). The per-dataset plots are left out; to include
# them again, add their sets here as
#   (folder in the Overleaf project, source folder of output/plots, mode text), e.g.
#   ("perEdge", OUT / "plots" / "perEdge", "perEdge"),
#   ("perBatch", OUT / "plots" / "perBatch", "perBatch"),
#   ("perBatchHPRC", OUT / "plots" / "perBatch" / "HPRC", "perBatch, HPRC")
FIGURE_SETS = []
FIGURE_WIDTH = 1600          # pixels: enough for the text width, keeps the project small
PLOT_TEXT = {
    "total_time_line": "Total time of every update (perEdge: an edge, perBatch: a haplotype) for every algorithm (log scale)",
    "total_time_bar": "Total time of every algorithm, summed over all updates",
    "components": "Time of every step (UpdateDS, DFS, TopoOrder, MinI, MaxO, SBCompute) over the updates, one panel per step (log scale)",
    "stacked_breakdown": "The first updates, every bar split into the time of the steps (right panel: algorithms below 1/20 of the largest, own scale)",
    "aggregated_versions": "Time of every step summed over all updates, one bar per algorithm (right panel: algorithms below 1/20 of the largest, own scale)",
    "activity": "Part of the graph touched by the updates: affected (reordered) and explored (reached by the searches) vertices and edges",
}
PLOT_ORDER = ["total_time_line", "total_time_bar", "components", "stacked_breakdown", "aggregated_versions", "activity"]


def shrink(src, dst):
    """copy a plot, at most FIGURE_WIDTH pixels wide"""
    dst.parent.mkdir(parents=True, exist_ok=True)
    with Image.open(src) as im:
        im = im.convert("RGB")
        if im.width > FIGURE_WIDTH:
            im = im.resize((FIGURE_WIDTH, round(im.height * FIGURE_WIDTH / im.width)), Image.LANCZOS)
        # 256-colour palette: plots use few colours, the files get several times smaller
        im = im.quantize(colors=256, method=Image.Quantize.MEDIANCUT)
        im.save(dst, optimize=True)


def dataset_order(name):
    """Ecoli first, then the chromosomes by number, X, Y"""
    if name.lower() == "ecoli":
        return (0, 0, "")
    m = re.match(r"chromosome_(\d+)$", name)
    return (1, int(m.group(1)), "") if m else (2, 0, name)


def build_figures():
    """copy (shrunk) the plots into figures/ and write appendix_figures.tex"""
    FIG = DEST / "figures"
    lines = ["% Appendix with the comparison plots (generated by Preprocessing/makeOverleaf.py), shrunk to",
             f"% {FIGURE_WIDTH} px width. Include after the tables:  \\input{{tables/appendix_figures}}",
             "% Needs \\usepackage{graphicx}. Move a plot into the text by cutting its figure block.", ""]
    count = 0

    def figure(rel, caption, label):
        nonlocal count
        lines.extend([r"\begin{figure}[ht]", r"\centering",
                      rf"\includegraphics[width=\linewidth]{{{rel}}}",
                      rf"\caption{{{caption}}}", rf"\label{{{label}}}", r"\end{figure}"])
        count += 1
        if count % 2 == 0:
            lines.append(r"\clearpage")

    for cls, cls_text in CLASSES:
      comp = sorted((OUT / "plots" / "comparison" / cls).glob("*.png"))
      hist = sorted((OUT / "plots" / "history" / cls).glob("*.png"))
      if comp or hist:
        lines += [rf"\section{{Plots: {cls_text}}}", ""]
      for f in hist:
        kind = next((k for k in HISTORY_TEXT if f.stem.endswith(k)), None)
        if kind is None:
            continue
        shrink(f, FIG / "history" / cls / f.name)
        figure(f"figures/history/{cls}/{f.name}", f"{HISTORY_TEXT[kind]} ({cls_text}).",
               f"fig:history:{cls}:{f.stem}")
      if comp:
        for f in comp:
            shrink(f, FIG / "comparison" / cls / f.name)
            mode = "perEdge" if "_perEdge_" in f.name else "perBatch"
            sets = f.stem.split(f"_{mode}_", 1)[-1]
            # a tag at the end (plotChromosomeComparison.py --tag): the data
            # structures for minI / maxO get their own caption
            subject = "one bar per algorithm"
            if sets.endswith("_scaling") and mode == "perBatch":
                figure(f"figures/comparison/{cls}/{f.name}",
                       "Total time against the number of inserted edges over all haplotypes, HPRC release 1 "
                       "chromosomes and release 2 chr22 (chr22 v2); all points measured. Left: our algorithms "
                       "(static, BRec, Batch, BatchH1) on a linear axis; right: BubbleGun and BubbleFinder on a "
                       "logarithmic axis, with BatchH1 as the reference. Our algorithms: sum of the timed "
                       "steps; the tools: reading the graph and finding the bubbles after every haplotype. "
                       "The time limit (dashed line) is not reached by any of them.",
                       f"fig:comparison:{cls}:{f.stem}")
            if sets.endswith("_scaling"):
                figure(f"figures/comparison/{cls}/{f.name}",
                       "Total time of the per-edge algorithms (log scale, sum of the timed steps) against the "
                       "number of inserted edges, HPRC release 1 chromosomes and release 2 chr22 (chr22 v2). "
                       "Filled: measured; hollow: static on release 2 as a lower bound: its measured time on "
                       "release 1 chr22 scaled by the graph size summed over the edge updates (static processes "
                       "the whole graph after every edge). No point: did not finish within the time limit "
                       "(dashed line). Only H1 finishes on every dataset. H1 spends 91--94\\% of its time on the "
                       "full recomputations of the superbubbles; an edge that goes against the order (case 4) "
                       "almost always forces one (98.7\\% on release 1 chr22), and only 5.9\\% of the edges of "
                       "release 2 are of this kind, against 17--21\\% on release 1. That is why release 2 chr22 "
                       "(3.3 million edges) takes less time than release 1 chr11 (2.7 million) and chrX "
                       "(2.4 million).", f"fig:comparison:{cls}:{f.stem}")
                continue
            if sets.endswith("_finished"):
                sets = sets[:-len("_finished")]
                subject = ("one bar per algorithm that finished on this dataset within the time limit (the "
                           "others are DNF, see the tables)")
            elif sets.endswith("_algorithms"):
                sets = sets[:-len("_algorithms")]
                subject = ("one bar per algorithm: static, the recomputation of minI / maxO (Rec), the "
                           "modified tree for minI (MinI) and the superbubble heuristic (H1)")
            elif sets.endswith("_trees"):
                sets = sets[:-len("_trees")]
                subject = ("one bar per data structure for minI and maxO, with all steps: static, normal "
                           "balanced binary search trees keyed by position (NBST), modified trees comparing "
                           "by the current positions (BST), modified tree for minI and maxO by comparison (MinI)")
            elif sets.endswith("_datastructures"):
                sets = sets[:-len("_datastructures")]
                subject = ("one bar per data structure for minI and maxO, only the time of maintaining "
                           "minI and maxO (the other steps are the same in all three: the same incremental "
                           "DFS and superbubble computation). NBST = normal balanced binary search trees "
                           "keyed by position, BST = modified trees comparing by the current positions "
                           "(minI and maxO), MinI = modified tree for minI and maxO by comparison. Zoom "
                           "panel: BST and MinI on their own scale")
            m = re.match(r"(HPRC|small|mixed)_(\d+)datasets$", sets)
            if m:   # short file name: describe the group
                sets = {"HPRC": f"all {m.group(2)} HPRC chromosomes",
                        "small": f"the {m.group(2)} small datasets (E. coli and human chromosomes)",
                        "mixed": f"{m.group(2)} datasets"}[m.group(1)]
            else:
                sets = sets.replace("_", ", ").replace("HPRCchr", "HPRC chr").replace("HPRCv2chr", "HPRC v2 chr")
            what = "Time of minI and maxO" if f.stem.endswith("_datastructures") else "Time of every step"
            dnf = (" A bar marked DNF: the run was not started because it was estimated to exceed the "
                   "time limit." if cls == "medium" and mode == "perEdge" else "")
            figure(f"figures/comparison/{cls}/{f.name}",
                   f"{what} summed over all updates ({mode}), one group of bars per dataset "
                   f"({sets}), {subject}.{dnf}", f"fig:comparison:{cls}:{f.stem}")
      if comp or hist:
        lines.append(r"\clearpage")

    for key, src, mode_text in FIGURE_SETS:
        if not src.is_dir():
            continue
        dsets = sorted([d for d in src.iterdir() if d.is_dir() and d.name not in ("HPRC", "HPRC2")], key=lambda d: dataset_order(d.name))
        if not dsets:
            continue
        lines += [rf"\section{{Plots ({mode_text})}}", ""]
        for d in dsets:
            label_ds = ("HPRC " if key == "perBatchHPRC" else "") + (d.name if d.name == "Ecoli" else d.name.replace("chromosome_", "chr"))
            lines += [rf"\subsection*{{{label_ds}}}"]
            plots = {}
            for f in d.glob("*.png"):
                kind = next((k for k in PLOT_ORDER if f.stem.endswith(k)), None)
                if kind:
                    plots[kind] = f
            for kind in PLOT_ORDER:
                if kind not in plots:
                    continue
                f = plots[kind]
                shrink(f, FIG / key / d.name / f.name)
                figure(f"figures/{key}/{d.name}/{f.name}", f"{label_ds} ({mode_text}): {PLOT_TEXT[kind]}.",
                       f"fig:{key}:{label_ds.replace(' ', '')}:{kind}")
            lines.append(r"\clearpage")
    (TABLES / "appendix_figures.tex").write_text("\n".join(lines) + "\n")
    return count


# what each paper table shows (for the catalogue)
PAPER_TEXT = {
    "table_component_summary": "time of every step summed over all datasets, speedup over static",
    "table_dataset_summary": "total time per dataset, speedup over static",
    "table_vertex_processed": "vertices explored / affected vs what static processes",
    "table_edge_processed": "edges explored / affected vs what static processes",
    "table_tool_bubblegun": "BubbleGun vs our algorithms (wall time, speedup = tool / algorithm)",
    "table_tool_bubblefinder": "BubbleFinder vs our algorithms (wall time, speedup = tool / algorithm)",
    "table_datastructures": "minI / maxO data structures: NBST vs BST vs MinI",
    "table_tools": "BubbleGun as reference: BubbleFinder, static and the best algorithm, speedup over BubbleGun",
    "table_full_time": "time of every step per dataset (one block per step), speedup over static",
    "table_full_time_edge": "long time table, perEdge: small / medium / large datasets, every step, speedup over static",
    **{f"table_{slug}_{m}": f"{name} time, {'perEdge' if m == 'edge' else 'perBatch'}: small / medium / large "
       f"datasets, static against the final algorithm" + (" (and the BST variant)" if slug == "maxo" else "")
       for slug, name in (("insert", "Insert"), ("dfs", "DFS"), ("toporder", "Top. Order"), ("mini", "MinI"),
                          ("maxo", "MaxO"), ("sbcompute", "SBCompute"), ("total", "Total")) for m in ("edge", "batch")},
    "table_full_time_batch": "long time table, perBatch: small / medium / large datasets, every step, speedup over static",
    "table_full_time_trees": "the same for the minI / maxO data structures: static, NBST, BST, MinI",
    "table_size": "datasets: haplotypes, base and final graph size",
    "table_dag_vs_original": "BubbleFinder on the original graphs and on the DAGs: superbubbles in both, only in each, touched sites",
    "table_datasets": "dataset table: haplotypes, total nodes and edges, small / medium / large",
    "table_memory_edge": "peak memory, perEdge: small / medium / large datasets, every algorithm, ratio to static",
    "table_memory_batch": "peak memory, perBatch: small / medium / large datasets, every algorithm, ratio to static",
    "table_tools_edge": "tools, perEdge: BubbleGun as reference, BubbleFinder, static and our algorithms, speedup over BubbleGun",
    "table_tools_batch": "tools, perBatch: BubbleGun as reference, BubbleFinder, static and our algorithms, speedup over BubbleGun",
    "table_step_UpdateDS": "time of UpdateDS per dataset, speedup over static",
    "table_step_DFS": "time of the DFS per dataset, speedup over static",
    "table_step_TopoOrder": "time of the topological order per dataset, speedup over static",
    "table_step_MinI": "time of minI per dataset, speedup over static",
    "table_step_MaxO": "time of maxO per dataset, speedup over static",
    "table_step_SBCompute": "time of the superbubble computation per dataset, speedup over static",
    "table_history_summary": "superbubble history per haplotype: base, born, promoted, died, final, DAG-affected",
    "table_history_var_types": "variant types of the final superbubbles",
    "table_edge_history": "superbubble history per edge, by case of the incremental DFS",
}


def convert(src, dst, prefix=None, tall=False):
    """copy a table: adjustbox around the tabular, label/ref prefix (detailed tables)"""
    s = src.read_text()
    box = r"\begin{adjustbox}{max width=\linewidth" + (r", max totalheight=0.85\textheight" if tall else "") + "}"
    if r"\begin{adjustbox}" not in s and r"\resizebox" not in s:
        s = s.replace(r"\begin{tabular}", box + "\n" + r"\begin{tabular}", 1)
        s = s.replace(r"\end{tabular}", r"\end{tabular}" + "\n" + r"\end{adjustbox}", 1)
    if prefix:
        s = re.sub(r"\\(label|ref)\{tab:", lambda m: f"\\{m.group(1)}{{tab:{prefix}:", s)
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(s)


# tables not copied; they stay in output/paperTables and their captions in
# table_captions.md: the per-step tables (each a block of table_full_time).
# The two total-time tables (per dataset: table_dataset_summary; per step:
# table_component_summary) are always kept. Where most runs did not finish
# (medium and large perEdge) the long table is left out too; the long table of
# medium perEdge is the one of chr22 alone (mediumPerEdgeChr22, of which only
# that table is copied).
LEFT_OUT = ("table_step_",)
DNF_MODES = ("mediumPerEdge", "largePerEdge")


def is_left_out(mode, stem):
    if stem == "table_full_time":
        return True                                     # replaced by the grouped long tables
    if mode != "grouped" and stem in ("table_vertex_processed", "table_edge_processed", "table_tools",
                                      "table_memory"):
        return True                                     # replaced by the grouped tables (perEdge and perBatch)
    if mode == "mediumPerEdgeChr22":
        return stem != "table_full_time_trees"
    if mode in DNF_MODES:
        return stem.startswith(("table_step_", "table_full_time"))
    return stem.startswith(LEFT_OUT)


def caption_of(path):
    """the caption text of a table file (the braces of \\caption{...} balanced)"""
    s = path.read_text()
    i = s.find(r"\caption{")
    if i < 0:
        return ""
    i += len(r"\caption{")
    depth, j = 1, i
    while j < len(s) and depth:
        depth += {"{": 1, "}": -1}.get(s[j], 0)
        j += 1
    return s[i:j - 1].strip()


def write_captions(captions):
    """table_captions.md: every paper table of every class and mode with its
    file, label, whether it is in the Overleaf folder, and its caption"""
    lines = ["# Captions of all paper tables", "",
             "Generated by Preprocessing/makeOverleaf.py. Every table of output/paperTables, also the ones "
             "left out of the Overleaf folder: the per-step tables (each a block of table_full_time); in "
             "medium and large "
             "perEdge, where most runs did not finish, the long table instead (medium: the long table of chr22 "
             "alone is in). A left-out table is still in output/paperTables/<class>/<mode>/.",
             ""]
    for mode in [m for m, _ in PAPER]:
        rows = [c for c in captions if c[0] == mode]
        if not rows:
            continue
        lines += [f"## {MODE_TITLE[mode]}", ""]
        for _, f, left_out in rows:
            lines += [f"### {f.stem}" + ("  (left out of Overleaf)" if left_out else ""), "",
                      f"- file: `{f.relative_to(REPO)}`", f"- label: `{label_of(f)}`",
                      f"- what: {PAPER_TEXT.get(f.stem, f.stem)}", "", caption_of(f), ""]
    text = "\n".join(lines) + "\n"
    (DEST / "table_captions.md").write_text(text)
    (OUT / "paperTables" / "CAPTIONS.md").write_text(text)


def label_of(path):
    m = re.search(r"\\label\{([^}]*)\}", path.read_text())
    return m.group(1) if m else ""


def main():
    if DEST.exists():
        shutil.rmtree(DEST)
    TABLES.mkdir(parents=True)

    # detailed tables (appendix)
    present = {}
    for mode, src in DETAILED:
        if not src.is_dir():
            continue
        for f in sorted(src.glob("*.tex")):
            if mode == "perBatch" and f.stem == "table_size_haplotypes":
                continue                                  # the same as perEdge
            convert(f, TABLES / mode / f.name, prefix=mode, tall=f.stem == "table_full_time")
            present.setdefault(mode, set()).add(f.stem)

    # paper tables (without the ones contained in table_full_time); every
    # caption, also of the left-out tables, goes into table_captions.md
    paper, captions = [], []
    for mode, src in PAPER:
        if not src.is_dir():
            continue
        for f in sorted(src.glob("*.tex")):
            left_out = is_left_out(mode, f.stem)
            captions.append((mode, f, left_out))
            if left_out:
                continue
            dst = TABLES / "paper" / mode / f.name
            convert(f, dst)
            paper.append((mode, f.stem, label_of(dst)))
    write_captions(captions)

    # tool comparison tables
    tools = []
    for t in TOOL_DIRS:
        for name in ("performance_perEdge", "performance_perBatch", "agreement"):
            f = OUT / "toolComparison" / t / f"{name}.tex"
            if f.exists():
                convert(f, TABLES / "toolComparison" / t / f.name, prefix="toolComparison-" + t.replace("_", "-"))
                tools.append((t, name))

    # catalogue of the paper tables
    lines = ["% Paper tables (generated by Preprocessing/makeOverleaf.py).",
             "% Copy each \\input line to the place in the text where the table belongs;",
             "% refer to it with \\ref{<label>}. Needs \\usepackage{adjustbox} in the preamble.", ""]
    for mode in [m for m, _ in PAPER]:
        rows = [p for p in paper if p[0] == mode]
        if not rows:
            continue
        lines.append(f"% ---- {MODE_TITLE[mode]} ----")
        for _, stem, label in rows:
            lines.append(f"% {PAPER_TEXT.get(stem, stem)}   (label: {label})")
            lines.append(f"\\input{{tables/paper/{mode}/{stem}}}")
        lines.append("")
    (TABLES / "paper_tables.tex").write_text("\n".join(lines) + "\n")

    # appendix
    order = ["table_size_haplotypes", "table_cases", "table_activity", "table_activity_totals_vertices",
             "table_activity_totals_edges", "table_batch", "table_TotalTime", "table_PeakMemory",
             "table_full_time", "table_UpdateDS", "table_DFS", "table_TopoOrder", "table_MinI",
             "table_MaxO", "table_SBCompute", "table_case_time", "table_case_time_best"]
    a = ["% Appendix with ALL tables (generated by Preprocessing/makeOverleaf.py): first the",
         "% paper tables (move them one by one into the text: cut the \\input line here and",
         "% paste it there), then the detailed tables.",
         "% Include at the end of the paper:  \\clearpage \\appendix \\input{tables/appendix_tables}",
         "% Needs \\usepackage{adjustbox} in the preamble.", ""]
    for mode in [m for m, _ in PAPER]:
        rows = [p for p in paper if p[0] == mode]
        if not rows:
            continue
        a += [f"\\section{{Paper tables ({MODE_TITLE[mode]})}}", f"\\label{{app:paper_{mode.lower()}}}"]
        for i, (_, stem, label) in enumerate(rows):
            a.append(f"% {PAPER_TEXT.get(stem, stem)}   (label: {label})")
            a.append(f"\\input{{tables/paper/{mode}/{stem}}}")
            if i % 2 == 1:
                a.append("\\clearpage")
        a += ["\\clearpage", ""]
    for mode in [m for m, _ in DETAILED]:
        stems = present.get(mode, set())
        if not stems:
            continue
        a += [f"\\section{{Detailed results ({MODE_TITLE[mode]})}}", f"\\label{{app:results_{mode.lower()}}}"]
        if mode == "perBatch":
            a.append("% Datasets: the same as in perEdge (Table~\\ref{tab:perEdge:size_haplotypes}).")
        rest = sorted(s for s in stems if s not in order)
        for i, stem in enumerate([s for s in order if s in stems] + rest):
            a.append(f"\\input{{tables/{mode}/{stem}}}")
            if i % 2 == 1 or stem == "table_full_time":
                a.append("\\clearpage")
        a += ["\\clearpage", ""]
    if tools:
        a += ["\\section{Comparison with the external tools}", "\\label{app:tools}"]
        for t, name in tools:
            a.append(f"\\input{{tables/toolComparison/{t}/{name}}}")
        a += ["\\clearpage", ""]
    (TABLES / "appendix_tables.tex").write_text("\n".join(a) + "\n")

    # plots
    n_fig = build_figures()
    print(f"figures: {n_fig} (appendix_figures.tex)")

    # zip
    zip_base = OUT / "overleaf_tables"
    if (OUT / "overleaf_tables.zip").exists():
        (OUT / "overleaf_tables.zip").unlink()
    shutil.make_archive(str(zip_base), "zip", DEST)
    n = sum(1 for _ in TABLES.rglob("*.tex"))
    print(f"{n} files in {TABLES}; zip: {zip_base}.zip")
    print(f"paper tables: {len(paper)}; detailed: {sum(len(v) for v in present.values())}; tool tables: {len(tools)}")


if __name__ == "__main__":
    main()
