#!/usr/bin/env python3
"""
edgeHistory.py - the per-edge superbubble history (what every inserted edge
destroys and creates) with the reference position and the haplotype support
of every edge, a summary per case, and figure 4.

The history itself comes from the C++ program (perEdge mode, third argument
"history"; see lib/findSuperbubble/superbubbleEdgeHistory.h):

    ./test/incMinIBSTH1 <dir>/incMinIBSTH1_time.txt perEdge history < numeric_edgelist.txt

=====================================================================
HOW TO RUN
=====================================================================

Run from the repository root:

    python3 Preprocessing/edgeHistory.py \\
        --events  output/HPRC/chromosome_22/history/incMinIBSTH1_edge_events.tsv \\
        --nodemap input/HPRC/chromosome_22/Graphs/chromosome_22.nodemap.tsv \\
        --graph   input/HPRC/chromosome_22/Graphs/numeric_edgelist.txt \\
        --gfa     input/HPRC/chromosome_22/Graphs/chromosome_22.gfa[.gz] \\
        --final   output/HPRC/chromosome_22/history/incMinIBSTH1_superbubbles.txt \\
        --label   "HPRC chr22 (release 1)" \\
        --out     output/HPRC/chromosome_22/history/chr22

Options:
    --events   <ID>_edge_events.tsv of the C++ program (required)
    --nodemap  <name>.nodemap.tsv of mcToDag.py (reference positions)
    --graph    numeric_edgelist.txt (the final DAG: a DFS order for the
               reference position of vertices off the reference)
    --gfa      the original GFA (haplotype walks: edge support)
    --final    <ID>_superbubbles.txt of the same run: check that replaying
               the edges backwards from its superbubbles with an inside
               vertex is consistent (needs --graph) (only for a finished
               run; without it, a run still going or stopped by its time limit
               is analysed up to its last complete edge)
    --top      edges listed in <out>_top_edges.tsv (default: 100)
    --label    title of figure 4
    --out      output prefix (required)

=====================================================================
OUTPUT
=====================================================================

Only superbubbles with at least one inside vertex are counted (trivial
superbubbles, a single edge with nothing inside, are not variant sites; see
lib/findSuperbubble/superbubbleEdgeHistory.h).

<out>_edges.tsv, one row per inserted edge that destroyed or created a
superbubble:
    edge_idx, step, genome, u, v, case, full, n_destroyed, n_created
                   as in the events
    ref_pos        position of u on the reference (0-based); for a vertex
                   off the reference the end of the nearest reference vertex
                   before it in a DFS order of the final graph; NA without
                   --nodemap / --graph
    ref_exact      1 if u itself lies on the reference
    edge_support   haplotypes whose walk uses the edge (with --gfa; NA
                   without)
<out>_top_edges.tsv   the --top edges with the most destroyed superbubbles
<out>_edge_summary.tsv, one row per case (and "all"): edges, edges that
    change something, destroyed, created, full recomputations
<out>_edge_destroyed.png   figure 4: superbubbles destroyed per edge along
    the reference, coloured by the edge's haplotype support
"""
import argparse
import collections
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mcToDag import opener, path_steps                                 # noqa: E402
from superbubbleHistory import read_nodemap, read_edges, dfs_order     # noqa: E402


def read_events(path):
    """rows of the edge events as dicts (lists left as strings); a last line
    without its newline (a run still writing, or stopped) is left out"""
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if not line.endswith("\n") or len(p) != len(head):
                print(f"[!] incomplete line left out (edge {p[0]}): the run is still writing or was stopped")
                continue
            yield dict(zip(head, p))


def check_replay(events, final_path, pos):
    """replays the edges backwards from the superbubbles of the final list
    that have an inside vertex (range entrance .. exit longer than the two in
    the DFS order pos of the final graph): every created superbubble must be
    alive and every destroyed one not; returns (inconsistencies,
    superbubbles of the base graph, trivial ones left out of the list)"""
    cur, trivial = set(), 0
    with open(final_path) as f:
        for line in f:
            p = line.split()
            if len(p) == 2 and p[0].isdigit() and (p[0], p[1]) != ("0", "1"):
                if pos[int(p[1])] - pos[int(p[0])] > 1:
                    cur.add(f"{p[0]}-{p[1]}")
                else:
                    trivial += 1
    bad = 0
    for r in reversed(events):
        if int(r["alive"]) != len(cur):
            bad += 1
        for c in ([] if r["created"] == "-" else r["created"].split(",")):
            bad += c not in cur
            cur.discard(c)
        for d in ([] if r["destroyed"] == "-" else r["destroyed"].split(",")):
            bad += d in cur
            cur.add(d)
    return bad, len(cur), trivial


def reference_positions(nodes, graph_path):
    """edge-list id -> (reference position, exact): on the reference its own
    position, else the end of the nearest reference vertex before it in a
    DFS order of the final graph"""
    pos, order = dfs_order(read_edges(graph_path))
    ref = {}
    last = None
    for x in order:
        n = nodes.get(x - 2)                       # edge list id = new id + 2
        if n is not None and n[5] and n[4] is not None:
            ref[x] = (n[4], 1)
            last = n[4] + n[2]
        elif last is not None:
            ref[x] = (last, 0)
    return ref


def edge_support(gfa, wanted, nodemap_path):
    """edge (u, v) in edge-list ids -> number of haplotypes whose walk uses it
    (in either direction of the walk), for the edges in 'wanted'"""
    node = {}
    with open(nodemap_path) as f:
        head = f.readline().rstrip("\n").split("\t")
        c = {n: i for i, n in enumerate(head)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            new = int(p[c["new"]])
            if new >= 0:
                node[p[c["orig"]]] = new + 2
    haps = collections.defaultdict(int)            # edge -> bit mask of haplotypes
    hap_ids = {}
    with opener(gfa) as fh:
        for line in fh:
            if line[0] not in "WP":
                continue
            _, hap, steps = path_steps(line.rstrip("\n").split("\t"))
            bit = 1 << hap_ids.setdefault(hap, len(hap_ids))
            w = [node.get(nm, -1) for nm, _ in steps]
            for a, b in zip(w, w[1:]):
                if a < 0 or b < 0:
                    continue
                e = (a, b) if (a, b) in wanted else (b, a) if (b, a) in wanted else None
                if e is not None:
                    haps[e] |= bit
    print(f"haplotypes in the walks: {len(hap_ids)}")
    return {e: bin(m).count("1") for e, m in haps.items()}


def main():
    ap = argparse.ArgumentParser(description="Per-edge superbubble history: tables and figure 4")
    ap.add_argument("--events", required=True)
    ap.add_argument("--nodemap")
    ap.add_argument("--graph")
    ap.add_argument("--gfa")
    ap.add_argument("--final")
    ap.add_argument("--top", type=int, default=100)
    ap.add_argument("--label", default="")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    if os.path.dirname(args.out):
        os.makedirs(os.path.dirname(args.out), exist_ok=True)

    events = list(read_events(args.events))
    print(f"edges: {len(events)}" + (f" (up to genome {events[-1]['step']}, {events[-1]['genome']})" if events else ""))
    pos = dfs_order(read_edges(args.graph))[0] if args.graph else None
    if args.final and pos is None:
        print("[!] the replay check needs --graph (which superbubbles of the final list have an inside vertex)")
    elif args.final:
        bad, n_base, n_triv = check_replay(events, args.final, pos)
        print(f"check: replay from {args.final} ({n_triv} trivial left out): {bad} inconsistencies, "
              f"{n_base} superbubbles with an inside vertex in the base graph")

    # ---- summary per case
    summ = collections.OrderedDict()
    for r in events:
        for k in (r["case"], "all"):
            s = summ.setdefault(k, [0, 0, 0, 0, 0])
            nd, nc = int(r["n_destroyed"]), int(r["n_created"])
            s[0] += 1
            s[1] += (nd + nc) > 0
            s[2] += nd
            s[3] += nc
            s[4] += r["full"] == "1"
    with open(args.out + "_edge_summary.tsv", "w") as o:
        o.write("case\tedges\tchanging\tdestroyed\tcreated\tfull\n")
        for k in sorted(summ, key=lambda k: (k == "all", k)):
            o.write(k + "\t" + "\t".join(map(str, summ[k])) + "\n")
    a = summ["all"]
    print(f"edges that change the superbubbles: {a[1]} of {a[0]} ({100.0 * a[1] / max(a[0], 1):.1f}%); "
          f"destroyed {a[2]}, created {a[3]}, full recomputations {a[4]}")

    # ---- edges that change something: reference position and support
    changing = [r for r in events if int(r["n_destroyed"]) + int(r["n_created"]) > 0]
    ref = {}
    if args.nodemap and args.graph:
        nodes, has_ref = read_nodemap(args.nodemap)
        if has_ref:
            ref = reference_positions(nodes, args.graph)
    support = {}
    if args.gfa and args.nodemap:
        support = edge_support(args.gfa, {(int(r["u"]), int(r["v"])) for r in changing}, args.nodemap)
    cols = ["edge_idx", "step", "genome", "u", "v", "case", "full", "n_destroyed", "n_created"]
    rows = []
    for r in changing:
        u, v = int(r["u"]), int(r["v"])
        rp, ex = ref.get(u, ("NA", "NA"))
        sp = support.get((u, v), 0) if support else "NA"
        rows.append([r[c] for c in cols] + [rp, ex, sp])
    head = "\t".join(cols + ["ref_pos", "ref_exact", "edge_support"]) + "\n"
    with open(args.out + "_edges.tsv", "w") as o:
        o.write(head)
        for row in rows:
            o.write("\t".join(map(str, row)) + "\n")
    top = sorted(rows, key=lambda row: -int(row[7]))[:args.top]
    with open(args.out + "_top_edges.tsv", "w") as o:
        o.write(head)
        for row in top:
            o.write("\t".join(map(str, row)) + "\n")
    print(f"written: {args.out}_edges.tsv ({len(rows)} edges), {args.out}_top_edges.tsv, {args.out}_edge_summary.tsv")

    # ---- figure 4: destroyed superbubbles per edge along the reference
    pts = [(row[9], int(row[7]), row[11]) for row in rows
           if row[9] != "NA" and int(row[7]) > 0]
    if not pts:
        print("[!] figure 4 needs --nodemap and --graph (reference positions)")
        return
    pts.sort(key=lambda p: p[2] if p[2] != "NA" else 0)      # well-supported edges drawn on top
    fig, ax = plt.subplots(figsize=(14, 4.5))
    x = [p[0] / 1e6 for p in pts]
    y = [p[1] for p in pts]
    if support:
        c = [max(p[2], 1) for p in pts]
        sc = ax.scatter(x, y, c=c, s=6, cmap="viridis", norm=mcolors.LogNorm(vmin=1, vmax=max(c)), linewidths=0)
        fig.colorbar(sc, ax=ax, pad=0.01).set_label("haplotypes using the edge")
    else:
        ax.scatter(x, y, s=6, color="#0072B2", linewidths=0)
    ax.set_yscale("log")
    ax.set_xlabel("position on the reference (Mb)")
    ax.set_ylabel("superbubbles destroyed by the edge")
    ax.set_title((args.label + ": " if args.label else "") +
                 f"superbubbles destroyed per inserted edge ({len(pts):,} edges destroy at least one)")
    ax.grid(True, axis="y", linestyle="--", alpha=0.4)
    fig.tight_layout()
    fig.savefig(args.out + "_edge_destroyed.png", dpi=200)
    print(f"written: {args.out}_edge_destroyed.png")


if __name__ == "__main__":
    main()
