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

<mode> = perEdge, perBatch (11 datasets), perBatchHPRC (HPRC chromosomes).
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

REPO = pathlib.Path(__file__).resolve().parent.parent
OUT = REPO / "output"
DEST = OUT / "overleaf_tables"
TABLES = DEST / "tables"

DETAILED = [("perEdge", OUT / "latexTablesperEdge"),
            ("perBatch", OUT / "latexTablesperBatch"),
            ("perBatchHPRC", OUT / "latexTablesperBatch_HPRC")]
PAPER = [("perEdge", OUT / "paperTables" / "perEdge"),
         ("perBatch", OUT / "paperTables" / "perBatch"),
         ("perBatchHPRC", OUT / "paperTables" / "perBatch_HPRC")]
TOOL_DIRS = ["bubblegun", "bubblefinder", "bubblegun_hprc22", "bubblefinder_hprc22",
             "bubblegun_hprc22_H1", "bubblefinder_hprc22_H1"]
MODE_TITLE = {"perEdge": "perEdge, 11 datasets", "perBatch": "perBatch, 11 datasets",
              "perBatchHPRC": "perBatch, HPRC chromosomes"}

# what each paper table shows (for the catalogue)
PAPER_TEXT = {
    "table_component_summary": "time of every step summed over all datasets, speedup over static",
    "table_dataset_summary": "total time per dataset, speedup over static",
    "table_vertex_processed": "vertices explored / affected vs what static processes",
    "table_edge_processed": "edges explored / affected vs what static processes",
    "table_tool_bubblegun": "BubbleGun vs our algorithms (wall time, speedup = tool / algorithm)",
    "table_tool_bubblefinder": "BubbleFinder vs our algorithms (wall time, speedup = tool / algorithm)",
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

    # paper tables
    paper = []
    for mode, src in PAPER:
        if not src.is_dir():
            continue
        for f in sorted(src.glob("*.tex")):
            dst = TABLES / "paper" / mode / f.name
            convert(f, dst)
            paper.append((mode, f.stem, label_of(dst)))

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
