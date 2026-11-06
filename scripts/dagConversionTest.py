"""
dagConversionTest.py

How much the conversion of an HPRC Minigraph-Cactus graph into a DAG
(Preprocessing/mcToDag.py) changes its superbubbles, for one chromosome of
HPRC release 1.1 (44 samples, 88 haplotypes, on GRCh38):

  superbubbles of the original graph (bidirected, BubbleFinder)
  against
  superbubbles of the DAG (directed, our algorithm ImpSB = incBSTBatchH2,
  after the last haplotype)

=====================================================================
HOW TO RUN
=====================================================================

Run from the repository root, with the virtual environment active. Needs
curl, vg (static binary from https://github.com/vgteam/vg/releases) and
BubbleFinder (see scripts/reviewerTest.py or the README).

    # chromosome 22, the smallest download
    python3 scripts/dagConversionTest.py -c 22 --vg ~/tools/vg

Choose the chromosome with -c. The five of the paper (download, share of the
original superbubbles that the DAG loses, measured for the paper):

    -c 22   0.94 GB   1.16 %
    -c 18   1.51 GB   0.00 %
    -c 17   1.54 GB   0.82 %
    -c X    1.83 GB   0.19 %
    -c 11   2.50 GB   0.20 %

Any other chromosome (1-22, X, Y) works the same way. The download is the
slow part: about 1 MB/s from the HPRC bucket (chr22: about 15 minutes); the
rest takes a few minutes. BubbleFinder needs up to about 7 GB of memory on
the original graph (chr11).

Options:
    -c / --chromosome   1 ... 22 X Y                      (required)
    --vg                vg program                        (default: vg)
    --threads           threads for vg convert and BubbleFinder (default: 4)
    --no-dag-check      do not run BubbleFinder on the DAG as well (the check
                        that our algorithm finds the same superbubbles there)
    --overwrite         run even if our results exist in output/HPRC/

=====================================================================
WHAT IT DOES
=====================================================================

  1. scripts/hprcPrepare.py --release 1 -c <c>: downloads chr<c>.vg, converts
     it to GFA (input/HPRC/chromosome_<c>/Graphs/chromosome_<c>.gfa.gz) and
     into a DAG with mcToDag.py (numeric_edgelist.txt, chromosome_<c>.dag.gfa,
     chromosome_<c>.nodemap.tsv = DAG vertex -> original segment,
     chromosome_<c>.stats.json = everything the conversion removed)
  2. BubbleFinder on the original GFA (bidirected)
  3. our algorithm on the DAG: scripts/superbubbleRunner.py -d HPRC
     -a incBSTBatchH2 -m perBatch (all haplotypes, one run); its superbubbles
     of the final graph
  4. (unless --no-dag-check) BubbleFinder on the DAG (directed): must give
     the same superbubbles as our algorithm
  5. tables: output/reviewer/dag_conversion_chr<c>.md (printed at the end)

A superbubble is compared as the unordered pair of its two ends, as original
segment names (the DAG's vertices are mapped back with the node map). The
superbubble at the artificial source and sink of our algorithm is left out.
"""
import argparse
import csv
import gzip
import json
import os
import pathlib
import shutil
import subprocess
import sys
import time

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "Preprocessing"))
from superbubbleTools.bubblefinderTool import DEFAULT_BINARY   # noqa: E402
from reviewerTest import log, run, table                       # noqa: E402

ALGO = "incBSTBatchH2"
OUT = REPO / "output" / "reviewer"


def bubblefinder(bf, gfa, out, directed, threads):
    """BubbleFinder's superbubbles: a set of unordered pairs of names as BubbleFinder writes them"""
    cmd = [bf, "superbubbles"] + (["--directed", "--gfa-directed"] if directed else ["--gfa"]) + \
          ["-j", str(threads), "-g", str(gfa), "-o", str(out)]
    log("$ " + " ".join(cmd))
    with open(str(out) + ".log", "w") as lf:
        if subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT).returncode != 0:
            sys.exit(f"[!] BubbleFinder failed, see {out}.log")
    with open(out) as f:
        n = int(f.readline())
        pairs = [tuple(line.split()) for line in f if line.strip()]
    assert len(pairs) == n, (out, len(pairs), n)
    return pairs


def main():
    ap = argparse.ArgumentParser(description="DAG conversion test on one HPRC release 1 chromosome")
    ap.add_argument("-c", "--chromosome", required=True,
                    help="1 ... 22, X or Y; the paper's: 22 (smallest), 18, 17, X, 11")
    ap.add_argument("--vg", default="vg")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--no-dag-check", action="store_true")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()
    c = args.chromosome.replace("chromosome_", "").replace("chr", "")
    name = f"chromosome_{c}"
    graphs = REPO / "input" / "HPRC" / name / "Graphs"
    ours_dir = REPO / "output" / "HPRC" / name / "timeFile" / "perbatch"
    work = OUT / f"dag_chr{c}"
    start = time.time()

    if (ours_dir / f"{ALGO}_time.txt").exists() and not args.overwrite:
        sys.exit(f"[!] {ours_dir.relative_to(REPO)} has results already; they would be replaced. "
                 f"Use --overwrite to go on.")
    bf = os.environ.get("BUBBLEFINDER", DEFAULT_BINARY)
    if not os.access(bf, os.X_OK):
        sys.exit(f"[!] BubbleFinder not found at {bf} (set BUBBLEFINDER=/path/to/BubbleFinder; install: README)")
    if shutil.which(args.vg) is None and not os.access(args.vg, os.X_OK):
        sys.exit(f"[!] vg not found ({args.vg}); give it with --vg (https://github.com/vgteam/vg/releases)")
    work.mkdir(parents=True, exist_ok=True)
    py = sys.executable

    log(f"1. download chr{c} (HPRC release 1.1) and convert it into a DAG")
    run([py, "-u", "scripts/hprcPrepare.py", "--release", "1", "-c", c, "--vg", args.vg,
         "--threads", str(args.threads)])

    log("2. BubbleFinder on the original graph")
    gfa = graphs / f"{name}.gfa"
    tmp = None
    if not gfa.exists():                                   # hprcPrepare.py keeps it gzipped
        tmp = work / "original.gfa"
        with gzip.open(graphs / f"{name}.gfa.gz", "rb") as fin, open(tmp, "wb") as fout:
            shutil.copyfileobj(fin, fout)
        gfa = tmp
    orig = bubblefinder(bf, gfa, work / "original_superbubbles.txt", False, args.threads)
    if tmp:
        tmp.unlink()

    log("3. our algorithm on the DAG")
    run([py, "scripts/superbubbleRunner.py", "-d", "HPRC", "-c", name, "-m", "perBatch", "-a", ALGO, "--runs", "1"])

    # DAG vertex -> original segment (node map of mcToDag.py; the edge list numbers vertex v as v + 2,
    # 0 and 1 are the artificial source and sink)
    new2orig = {}
    with open(graphs / f"{name}.nodemap.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if int(row["new"]) >= 0:
                new2orig[int(row["new"])] = row["orig"]
    ours = set()
    for line in open(ours_dir / f"{ALGO}_superbubbles.txt"):
        p = line.split()
        if len(p) < 2 or line.startswith("#"):
            continue
        a, b = int(p[0]), int(p[1])
        if a > 1 and b > 1:
            ours.add(frozenset((new2orig[a - 2], new2orig[b - 2])))
    O = {frozenset(p) for p in orig}

    check = None
    if not args.no_dag_check:
        log("4. BubbleFinder on the DAG (check)")
        dag = bubblefinder(bf, graphs / f"{name}.dag.gfa", work / "dag_superbubbles.txt", True, args.threads)
        D = {frozenset((new2orig[int(a.rstrip("+-"))], new2orig[int(b.rstrip("+-"))])) for a, b in dag}
        check = (D == ours, len(D))

    log("5. tables")
    s = json.load(open(graphs / f"{name}.stats.json"))

    def pct(sec, key):
        return f"{s[sec][key]:.3f} %"
    conv = table(["", f"chr{c}"], [
        ["Segments (original)", s["input"]["nodes"]], ["Links (original)", s["input"]["links"]],
        ["Haplotype walks", s["input"]["haplotype_walks"]],
        ["Segments flipped", f"{s['orientation']['nodes_flipped']} ({pct('orientation', 'pct_nodes_flipped')})"],
        ["Inversion links removed", f"{s['links_removed']['inversion_links']} ({pct('links_removed', 'pct_inversion')})"],
        ["Self-loops removed", s["links_removed"]["self_loops"]],
        ["Back edges removed (cycles)", s["links_removed"]["cycle_breaking_back_edges"]],
        ["Links removed with small components", s["links_removed"]["dropped_with_small_components"]],
        ["Links removed, total", f"{s['links_removed']['total_removed']} ({pct('links_removed', 'pct_total_removed')})"],
        ["Segments removed", f"{s['nodes_removed']['nodes']} ({pct('nodes_removed', 'pct_nodes')})"],
        ["Base pairs removed", f"{s['nodes_removed']['bp']} ({pct('nodes_removed', 'pct_bp')})"],
        ["Walk steps broken", f"{s['walks']['broken_transitions']} ({pct('walks', 'pct_broken')})"],
        ["DAG vertices", s["output"]["nodes"]], ["DAG edges", s["output"]["edges"]],
        ["DAG acyclic", "yes" if s["output"]["acyclic"] else "NO"]])
    sb = table(["", f"chr{c}"], [
        ["Superbubbles in the original (BubbleFinder)", len(O)],
        ["Superbubbles in the DAG (ImpSB)", len(ours)],
        ["In both", len(O & ours)],
        ["Only in the original (lost)", f"{len(O - ours)} ({100 * len(O - ours) / len(O):.2f} %)"],
        ["Only in the DAG", len(ours - O)],
        ["Preserved", f"{100 * len(O & ours) / len(O):.2f} %"]]
        + ([["BubbleFinder on the DAG gives the same superbubbles as ImpSB",
             ("yes" if check[0] else "NO") + f" ({check[1]})"]] if check else []))
    text = "\n\n".join([f"# DAG conversion of HPRC release 1.1, chr{c}",
                        "## What the conversion removed (mcToDag.py)", conv,
                        "## Superbubbles: original graph against the DAG "
                        "(unordered pairs of original segment names)", sb]) + "\n"
    out = OUT / f"dag_conversion_chr{c}.md"
    out.write_text(text)
    print("\n" + text)
    log(f"done in {(time.time() - start) / 60:.1f} min; tables in {out.relative_to(REPO)}")


if __name__ == "__main__":
    main()
