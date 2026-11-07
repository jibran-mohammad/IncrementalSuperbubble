#!/usr/bin/env python3
"""
mcToDag.py - Normalize a Minigraph-Cactus GFA (bidirected, possibly cyclic)
into a directed acyclic graph, report exactly what was removed, and
(optionally) write the input files of this repository's pipeline.

Experimental: for large pangenome graphs (e.g. HPRC Minigraph-Cactus
chromosomes), dataset "HPRC" in the runners (input/HPRC/chromosome_<X>/).

=====================================================================
HOW TO RUN
=====================================================================

    # 1. get a graph, e.g. HPRC chromosome 22 (.vg), and convert it to GFA
    #    (vg: static binary from https://github.com/vgteam/vg/releases)
    vg convert -f chr22.vg > input/HPRC/chromosome_22/Graphs/chromosome_22.gfa

    # 2. the reference sample name (first column after W)
    grep "^W" input/HPRC/chromosome_22/Graphs/chromosome_22.gfa | cut -f2 | sort -u | head

    # 3. normalize to a DAG and write the pipeline files next to it
    python3 Preprocessing/mcToDag.py input/HPRC/chromosome_22/Graphs/chromosome_22.gfa \
        --ref GRCh38 --pipeline-dir input/HPRC/chromosome_22/Graphs

    # (pipelineRunner.py -d HPRC runs this step itself; pypy3 is several
    #  times faster than python3 on large graphs)

Options:
    gfa               input GFA (.gfa or .gfa.gz) with W or P lines
    --ref             reference sample name (W sample field or P name prefix)
    --out             prefix of the normalization outputs below
                      (default: <pipeline-dir>/<name> if --pipeline-dir is given)
    --pipeline-dir    also write the pipeline input files into this folder
    --name            graph name in numeric_edgelist.txt (default: the GFA
                      file name without .gfa / .gfa.gz)
    --all-components  keep all weakly connected components (default: only
                      the one containing the start of the reference)

=====================================================================
WHAT IT DOES
=====================================================================

  1. Orient every node. Reference nodes take the orientation the reference
     walk uses; other nodes take the majority orientation over all haplotype
     walks (ties -> forward), each walk counted in the reference's
     direction: a walk that uses the reference nodes mostly in the opposite
     orientation (a reverse-complemented contig) votes reversed.
  2. Convert each bidirected link to a directed edge under that orientation.
     Links that would need one endpoint reversed (inversion links) cannot be
     represented in a directed graph and are dropped. Self-loops are dropped.
  3. Give each node a coordinate: its rank on the reference walk, or for
     non-reference nodes the average position interpolated between the
     flanking reference nodes in every walk that visits it.
  4. Break cycles with a DFS that starts at the reference start and explores
     successors in coordinate order; only DFS back edges are removed, so no
     edge that lies on no cycle is ever removed.
  5. Keep only the weakly connected component that contains the start of
     the reference (--all-components keeps all), verify the result is
     acyclic, and split the edges into batches: the reference's edges first,
     then for every haplotype (sample#hap, in file order) the NEW edges it
     contributes, in walk order. Edges that no walk uses are
     left out of the pipeline files (as in the rest of the pipeline).

=====================================================================
OUTPUT
=====================================================================

  Normalization (prefix --out):
    <out>.dag.gfa      S/L lines only, nodes renumbered 0..n-1, all links +/+
    <out>.batches.tsv  one batch per haplotype, reference first ("#" lines)
    <out>.removed_links.tsv  every removed link: kind (inversion, self_loop,
                       back_edge, component), its endpoints as original names
                       and as new ids (-1: the node was removed too); for
                       Preprocessing/superbubbleHistory.py --removed-links
                       (dag_affected)
    <out>.nodemap.tsv  original node name -> new id, flipped?, length, segment,
                       reference sequence and 0-based position (reference
                       nodes: their start; other nodes: the end of the
                       reference node before them), on the reference?
    <out>.stats.json   all counts and percentages

  Pipeline files (--pipeline-dir), the same as Preprocessing/ writes for the
  other datasets; node with new id k is segment s{k+1} and number k+2 in the
  edge list (0 and 1 are the algorithms' artificial source and sink):
    numeric_edgelist.txt          base graph = reference edges, then per
                                  haplotype its new edges (C++ algorithms,
                                  checkSuperbubbles.py)
    reduced_graph.gfa             base graph as GFA (BubbleGun)
    unique_elements.txt           per haplotype its new NODE / EDGE lines
    ordered_haplotypes.txt        order of the haplotypes (BubbleGun)
    haplotype_removal_report.txt  sizes and what was removed (summary of
                                  pipelineRunner.py)
"""
import argparse
import math
import gzip
import json
import re
import os
import sys
import time
from array import array

STEP_RE = re.compile(r'([><])([^><]+)')


def opener(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def path_steps(fields):
    """Return (sample, hap_key, steps) for a W or P line; steps = [(name, rev)]."""
    if fields[0] == 'W':
        sample, hap, seqid, walk = fields[1], fields[2], fields[3], fields[6]
        steps = [(m.group(2), m.group(1) == '<') for m in STEP_RE.finditer(walk)]
        return sample, f'{sample}#{hap}', steps
    name = fields[1]
    sample = name.split('#')[0]
    parts = name.split('#')
    hap_key = '#'.join(parts[:2]) if len(parts) >= 2 else name
    steps = [(tok[:-1], tok[-1] == '-') for tok in fields[2].split(',') if tok]
    return sample, hap_key, steps


def ref_path_origin(fields):
    """(sequence name, 0-based start) of a reference W or P line: W has them
    in its seqid / seq_start fields; a P name like GRCh38#0#chr22 or
    GRCh38#0#chr22[1000-2000] is split at its last '#' (start in brackets)."""
    if fields[0] == 'W':
        try:
            start = int(fields[4])
        except ValueError:
            start = 0
        return fields[3], start
    name = fields[1].split('#')[-1]
    start = 0
    if name.endswith(']') and '[' in name:
        name, rng = name[:-1].split('[', 1)
        try:
            start = int(rng.split('-')[0])
        except ValueError:
            start = 0
    return name, start


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def write_pipeline_files(args, stats, batch_order, batches, ref_haps, unwalked,
                         eu, ev, new_id, length, seconds):
    """
    The input files of the pipeline (see OUTPUT above), from the batches:
    base graph = the reference's edges, then every other haplotype with
    the new edges it contributes. Node with new id k: segment s{k+1},
    number k+2 in numeric_edgelist.txt.
    """
    d = args.pipeline_dir
    os.makedirs(d, exist_ok=True)
    seg = lambda x: f's{new_id[x] + 1}'
    num = lambda x: new_id[x] + 2

    base = [k for hap in batch_order if hap in ref_haps for k in batches[hap]]
    haps = [hap for hap in batch_order if hap not in ref_haps]
    base_nodes = []
    seen = set()
    for k in base:
        for x in (eu[k], ev[k]):
            if x not in seen:
                seen.add(x)
                base_nodes.append(x)
    n_base_nodes = len(base_nodes)

    # numeric_edgelist.txt (format of Preprocessing/edgeListProcessor.py)
    with open(os.path.join(d, 'numeric_edgelist.txt'), 'w') as o:
        o.write("# Node Numbering Convention: Node 'sX' is represented as integer 'X+1'.\n")
        o.write("# e.g., s1 -> 2, s2 -> 3, s10 -> 11.\n")
        o.write(f"# Written by Preprocessing/mcToDag.py from {args.gfa} (reference {args.ref}).\n\n")
        o.write(f'{args.name}\n{len(base)} {n_base_nodes}\n')
        for k in base:
            o.write(f'{num(eu[k])} {num(ev[k])}\n')
        for hap in haps:
            o.write(f'\n# HAPLOTYPE: {hap}\n')
            if batches[hap]:                     # no count line for a haplotype without new edges
                o.write(f'{len(batches[hap])}\n')
                for k in batches[hap]:
                    o.write(f'{num(eu[k])} {num(ev[k])}\n')

    # reduced_graph.gfa: the base graph
    with open(os.path.join(d, 'reduced_graph.gfa'), 'w') as o:
        o.write('H\tVN:Z:1.0\n')
        for x in base_nodes:
            o.write(f'S\t{seg(x)}\t*\tLN:i:{length[x]}\n')
        for k in base:
            o.write(f'L\t{seg(eu[k])}\t+\t{seg(ev[k])}\t+\t0M\n')

    # unique_elements.txt / ordered_haplotypes.txt (Preprocessing/haplotypeParser.py)
    with open(os.path.join(d, 'unique_elements.txt'), 'w') as o, \
         open(os.path.join(d, 'ordered_haplotypes.txt'), 'w') as oh:
        oh.write('# Haplotypes in the order they are added (Preprocessing/mcToDag.py)\n\n')
        for hap in haps:
            oh.write(f'# HAPLOTYPE: {hap}\n')
            o.write(f'# HAPLOTYPE: {hap}\n')
            for k in batches[hap]:
                for x in (eu[k], ev[k]):
                    if x not in seen:
                        seen.add(x)
                        o.write(f'NODE\t{seg(x)}\n')
            for k in batches[hap]:
                o.write(f'EDGE\t{seg(eu[k])}+->{seg(ev[k])}+\n')
            o.write('\n')

    # haplotype_removal_report.txt (keys read by pipelineRunner.py's summary)
    n_hap_edges = sum(len(batches[h]) for h in haps)
    lr, w, out = stats['links_removed'], stats['walks'], stats['output']
    with open(os.path.join(d, 'haplotype_removal_report.txt'), 'w') as o:
        o.write(f"graph: {os.path.abspath(args.gfa)}\n")
        o.write("prepared by: Preprocessing/mcToDag.py (Minigraph-Cactus GFA -> DAG)\n")
        o.write(f"walks (genomes, incl. reference): {len(haps) + 1}\n")
        o.write(f"reference walk: {args.ref} ({', '.join(sorted(ref_haps))})\n")
        o.write(f"segments in graph: {stats['input']['nodes']}\n")
        o.write(f"links in graph: {stats['input']['links']}\n")
        o.write(f"links removed (inversion links): {lr['inversion_links']}\n")
        o.write(f"links removed (self-loops): {lr['self_loops']}\n")
        o.write(f"links removed (cycle back edges): {lr['cycle_breaking_back_edges']}\n")
        o.write(f"links removed (small components): {lr['dropped_with_small_components']}\n")
        o.write(f"links used by no walk (dropped): {len(unwalked)}\n")
        o.write(f"walk transitions broken: {w['broken_transitions']} of {w['transitions']} ({w['pct_broken']}%)\n")
        o.write(f"acyclic: {'yes' if out['acyclic'] else 'NO'}\n")
        o.write(f"base graph segments (= reference walk): {n_base_nodes}\n")
        o.write(f"base graph links (= reference walk): {len(base)}\n")
        o.write(f"segments added by haplotypes: {len(seen) - n_base_nodes}\n")
        o.write(f"links added by haplotypes: {n_hap_edges}\n")
        o.write("element given to: first haplotype (file order) whose walk uses it\n")
        o.write(f"preprocessing time (s): {seconds:.3f}\n")
    log(f'pipeline files written to {d}')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('gfa')
    ap.add_argument('--ref', required=True,
                    help='reference sample name (W sample field or P name prefix), e.g. GRCh38')
    ap.add_argument('--out', help='output prefix (default: <pipeline-dir>/<name>)')
    ap.add_argument('--pipeline-dir', help='also write the pipeline input files into this folder')
    ap.add_argument('--name', help='graph name in numeric_edgelist.txt (default: GFA file name)')
    ap.add_argument('--all-components', action='store_true',
                    help='keep all weakly connected components (default: only the one containing the start of the reference)')
    args = ap.parse_args()
    if args.name is None:
        base_name = os.path.basename(args.gfa)
        for ext in ('.gz', '.gfa'):
            if base_name.endswith(ext):
                base_name = base_name[:-len(ext)]
        args.name = base_name
    if args.out is None:
        if args.pipeline_dir is None:
            ap.error('give --out or --pipeline-dir')
        args.out = os.path.join(args.pipeline_dir, args.name)
    for folder in (os.path.dirname(args.out), args.pipeline_dir):
        if folder:
            os.makedirs(folder, exist_ok=True)
    t_start = time.time()

    # ---------------- pass 1: nodes, links, reference walk ----------------
    log('pass 1: nodes, links, reference')
    ids = {}
    names = []
    length = array('q')

    def nid(name):
        i = ids.get(name)
        if i is None:
            i = len(names)
            ids[name] = i
            names.append(name)
            length.append(0)
        return i

    la, lb = array('i'), array('i')
    lra, lrb = bytearray(), bytearray()
    ref_rank = {}
    ref_rev = {}
    n_ref_paths = 0
    # reference coordinates: per rank the sequence name and the 0-based bp
    # position of the node on it (path start + lengths of the nodes before);
    # computed after pass 1, when all segment lengths are known
    ref_paths = []
    rank_seq, rank_pos, rank_node = [], array('q'), array('i')

    with opener(args.gfa) as fh:
        for line in fh:
            c = line[0]
            if c == 'S':
                f = line.rstrip('\n').split('\t')
                i = nid(f[1])
                if f[2] != '*':
                    length[i] = len(f[2])
                else:
                    for tag in f[3:]:
                        if tag.startswith('LN:i:'):
                            length[i] = int(tag[5:])
            elif c == 'L':
                f = line.split('\t', 5)
                la.append(nid(f[1])); lra.append(f[2] == '-')
                lb.append(nid(f[3])); lrb.append(f[4] == '-')
            elif c in 'WP':
                f = line.rstrip('\n').split('\t')
                sample, _, steps = path_steps(f)
                if sample != args.ref:
                    continue
                n_ref_paths += 1
                seq, start = ref_path_origin(f)
                path_nodes = array('i')
                for name, rev in steps:
                    i = nid(name)
                    path_nodes.append(i)
                    if i not in ref_rank:
                        ref_rank[i] = len(ref_rank)
                        ref_rev[i] = rev
                ref_paths.append((seq, start, path_nodes))
    n = len(names)
    m_in = len(la)
    rank_seq = [None] * len(ref_rank)
    rank_pos = array('q', bytes(8 * len(ref_rank)))
    rank_node = array('i', bytes(4 * len(ref_rank)))
    for seq, start, path_nodes in ref_paths:
        pos = start
        for i in path_nodes:
            r = ref_rank[i]
            if rank_seq[r] is None:          # first occurrence on the reference
                rank_seq[r], rank_pos[r], rank_node[r] = seq, pos, i
            pos += length[i]
    del ref_paths
    if not ref_rank:
        sys.exit(f'error: no W/P line with sample "{args.ref}" found; '
                 'check names with: grep "^W" file.gfa | cut -f2,3,4 | sort -u | head')
    log(f'  {n} nodes, {m_in} links, reference: {n_ref_paths} path(s), {len(ref_rank)} nodes')

    # ------------- pass 2: orientation votes + coordinates ----------------
    log('pass 2: haplotype walks (orientation, coordinates)')
    plus = array('i', bytes(4 * n))
    minus = array('i', bytes(4 * n))
    ksum = array('d', bytes(8 * n))
    kcnt = array('i', bytes(4 * n))
    n_walks = 0
    n_walks_reversed = n_walks_without_ref = 0
    with opener(args.gfa) as fh:
        for line in fh:
            if line[0] not in 'WP':
                continue
            f = line.rstrip('\n').split('\t')
            sample, _, steps = path_steps(f)
            if sample == args.ref:
                continue
            n_walks += 1
            st = [(ids[nm], rv) for nm, rv in steps if nm in ids]
            # Direction of the walk: many assembly contigs are stored reverse-
            # complemented and walk the graph backwards (<6<5<4 ...). On the
            # reference nodes such a walk mostly uses the opposite orientation
            # of the reference; its votes are then counted reversed, so that
            # every walk votes in the reference's direction. (Without this, a
            # node passed backwards by more reversed contigs than forward ones
            # was flipped, and its links were dropped as inversions.)
            same = opposite = 0
            for x, rv in st:
                if x in ref_rank:
                    if rv == ref_rev[x]:
                        same += 1
                    else:
                        opposite += 1
            walk_rev = opposite > same
            n_walks_reversed += walk_rev
            n_walks_without_ref += (same + opposite == 0)
            L = len(st)
            lc, li = [None] * L, [0] * L
            last_c, last_i = None, 0
            for i, (x, rv) in enumerate(st):
                if rv != walk_rev:
                    minus[x] += 1
                else:
                    plus[x] += 1
                r = ref_rank.get(x)
                if r is not None:
                    last_c, last_i = r, i
                lc[i], li[i] = last_c, last_i
            next_c, next_i = None, 0
            for i in range(L - 1, -1, -1):
                x = st[i][0]
                r = ref_rank.get(x)
                if r is not None:
                    next_c, next_i = r, i
                    continue
                a, ai = lc[i], li[i]
                if a is not None and next_c is not None and next_i != ai:
                    coord = a + (next_c - a) * (i - ai) / (next_i - ai)
                elif a is not None:
                    coord = a + 1e-3 * (i - ai)
                elif next_c is not None:
                    coord = next_c - 1e-3 * (next_i - i)
                else:
                    continue
                ksum[x] += coord
                kcnt[x] += 1

    flip = bytearray(n)
    key = array('d', bytes(8 * n))
    has_key = bytearray(n)
    for x in range(n):
        if x in ref_rank:
            flip[x] = ref_rev[x]
            key[x] = ref_rank[x]
            has_key[x] = 1
        else:
            flip[x] = minus[x] > plus[x]
            if kcnt[x]:
                key[x] = ksum[x] / kcnt[x]
                has_key[x] = 1

    # ---------------- orient links into directed edges --------------------
    log('orienting links')
    edge_idx = {}
    eu, ev = array('i'), array('i')
    n_inv = n_self = n_dup = 0
    dropped_links = []            # (kind, a, b) of the links that cannot be directed edges
    for k in range(m_in):
        a, b = la[k], lb[k]
        oa, ob = lra[k] ^ flip[a], lrb[k] ^ flip[b]
        if a == b:
            n_self += 1
            dropped_links.append(('self_loop', a, b))
            continue
        if not oa and not ob:
            u, v = a, b
        elif oa and ob:
            u, v = b, a
        else:
            n_inv += 1
            dropped_links.append(('inversion', a, b))
            continue
        if (u, v) in edge_idx:
            n_dup += 1
            continue
        edge_idx[(u, v)] = len(eu)
        eu.append(u); ev.append(v)
    m_dir = len(eu)

    # nodes never visited by any walk: borrow a coordinate from neighbours
    for _ in range(3):
        acc = {}
        for k in range(m_dir):
            u, v = eu[k], ev[k]
            if has_key[u] and not has_key[v]:
                s = acc.setdefault(v, [0.0, 0]); s[0] += key[u] + 0.5; s[1] += 1
            if has_key[v] and not has_key[u]:
                s = acc.setdefault(u, [0.0, 0]); s[0] += key[v] - 0.5; s[1] += 1
        if not acc:
            break
        for x, (s, c) in acc.items():
            key[x] = s / c
            has_key[x] = 1
    big = float(len(ref_rank) + 1)
    for x in range(n):
        if not has_key[x]:
            key[x] = big + x

    # ---------------- break cycles: DFS back edges ------------------------
    log('breaking cycles')
    adj = [[] for _ in range(n)]
    indeg = array('i', bytes(4 * n))
    for k in range(m_dir):
        adj[eu[k]].append(k)
        indeg[ev[k]] += 1
    for lst in adj:
        if len(lst) > 1:
            lst.sort(key=lambda k: key[ev[k]])

    removed = bytearray(m_dir)
    color = bytearray(n)
    ptr = array('i', bytes(4 * n))
    ref_start = min(ref_rank, key=ref_rank.get)
    roots = [ref_start] + sorted(range(n), key=lambda x: (indeg[x] > 0, key[x]))
    for r in roots:
        if color[r]:
            continue
        color[r] = 1
        stack = [r]
        while stack:
            x = stack[-1]
            p = ptr[x]
            if p < len(adj[x]):
                ptr[x] = p + 1
                k = adj[x][p]
                y = ev[k]
                if color[y] == 0:
                    color[y] = 1
                    stack.append(y)
                elif color[y] == 1:
                    removed[k] = 1
            else:
                color[x] = 2
                stack.pop()
    n_back = sum(removed)
    n_back_ref = sum(1 for k in range(m_dir) if removed[k]
                     and eu[k] in ref_rank and ev[k] in ref_rank
                     and ref_rank[ev[k]] == ref_rank[eu[k]] + 1)

    # ---------------- connected components --------------------------------
    log('components')
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for k in range(m_dir):
        if not removed[k]:
            a, b = find(eu[k]), find(ev[k])
            if a != b:
                parent[a] = b
    comp_size = {}
    for x in range(n):
        r = find(x)
        comp_size[r] = comp_size.get(r, 0) + 1
    n_comp = len(comp_size)
    if args.all_components:
        keep = bytearray([1]) * n
    else:
        best = find(ref_start)
        keep = bytearray(n)
        for x in range(n):
            keep[x] = find(x) == best
    for k in range(m_dir):
        if not removed[k] and not (keep[eu[k]] and keep[ev[k]]):
            removed[k] = 2  # dropped with its component

    # ---------------- renumber in coordinate order, verify DAG ------------
    kept_nodes = sorted((x for x in range(n) if keep[x]), key=lambda x: key[x])
    new_id = array('i', [-1]) * n
    for i, x in enumerate(kept_nodes):
        new_id[x] = i
    N = len(kept_nodes)
    kept_edges = [k for k in range(m_dir) if not removed[k]]
    out_adj = [[] for _ in range(N)]
    indeg2 = [0] * N
    for k in kept_edges:
        u, v = new_id[eu[k]], new_id[ev[k]]
        out_adj[u].append(v)
        indeg2[v] += 1
    q = [i for i in range(N) if indeg2[i] == 0]
    sources = len(q)
    seen = 0
    tmp = indeg2[:]
    while q:
        x = q.pop()
        seen += 1
        for y in out_adj[x]:
            tmp[y] -= 1
            if tmp[y] == 0:
                q.append(y)
    acyclic = seen == N
    sinks = sum(1 for i in range(N) if not out_adj[i])

    # ---------------- pass 3: per-haplotype batches ----------------------
    log('pass 3: per-haplotype edge batches')
    batch_order, batches = [], {}
    ref_haps = set()
    used = bytearray(m_dir)
    trans = broken = 0
    haps_broken = set()

    def add_walk(hap, steps):
        nonlocal trans, broken
        if hap not in batches:
            batches[hap] = []
            batch_order.append(hap)
        out = batches[hap]
        prev = None
        for nm, rv in steps:
            x = ids.get(nm)
            if x is None:
                prev = None
                continue
            ox = rv ^ flip[x]
            if prev is not None:
                px, po = prev
                trans += 1
                if px == x:
                    k = None
                elif not po and not ox:
                    k = edge_idx.get((px, x))
                elif po and ox:
                    k = edge_idx.get((x, px))
                else:
                    k = None
                if k is None or removed[k]:
                    broken += 1
                    haps_broken.add(hap)
                elif not used[k]:
                    used[k] = 1
                    out.append(k)
            prev = (x, ox)

    for want_ref in (True, False):
        with opener(args.gfa) as fh:
            for line in fh:
                if line[0] not in 'WP':
                    continue
                f = line.rstrip('\n').split('\t')
                sample, hap, steps = path_steps(f)
                if (sample == args.ref) == want_ref:
                    add_walk(hap, steps)
                    if want_ref:
                        ref_haps.add(hap)
    unwalked = [k for k in kept_edges if not used[k]]

    # ---------------- outputs ----------------------------------------------
    log('writing outputs')
    with open(f'{args.out}.dag.gfa', 'w') as o:
        o.write('H\tVN:Z:1.0\n')
        for x in kept_nodes:
            o.write(f'S\t{new_id[x]}\t*\tLN:i:{length[x]}\n')
        for k in kept_edges:
            o.write(f'L\t{new_id[eu[k]]}\t+\t{new_id[ev[k]]}\t+\t0M\n')
    with open(f'{args.out}.batches.tsv', 'w') as o:
        for hap in batch_order:
            o.write(f'#\t{hap}\t{len(batches[hap])}\n')
            for k in batches[hap]:
                o.write(f'{new_id[eu[k]]}\t{new_id[ev[k]]}\n')
        if unwalked:
            o.write(f'#\tunwalked\t{len(unwalked)}\n')
            for k in unwalked:
                o.write(f'{new_id[eu[k]]}\t{new_id[ev[k]]}\n')
    with open(f'{args.out}.removed_links.tsv', 'w') as o:
        o.write('kind\tu_orig\tv_orig\tu_new\tv_new\n')
        for kind, a, b in dropped_links:
            o.write(f'{kind}\t{names[a]}\t{names[b]}\t{new_id[a]}\t{new_id[b]}\n')
        for k in range(m_dir):
            if removed[k]:
                kind = 'back_edge' if removed[k] == 1 else 'component'
                o.write(f'{kind}\t{names[eu[k]]}\t{names[ev[k]]}\t{new_id[eu[k]]}\t{new_id[ev[k]]}\n')
    with open(f'{args.out}.nodemap.tsv', 'w') as o:
        o.write('orig\tnew\tflipped\tlength\tsegment\tref_seq\tref_pos\ton_ref\n')
        last_rank = len(rank_pos) - 1
        for x in range(n):
            seg = f's{new_id[x] + 1}' if new_id[x] >= 0 else '-'
            # reference coordinate: a reference node its own start; any other
            # node the end of the reference node before it; '-' if it has no
            # rank key
            r = ref_rank.get(x)
            if r is not None:
                seq, rpos, on_ref = rank_seq[r], rank_pos[r], 1
            elif has_key[x] and key[x] <= 0:
                seq, rpos, on_ref = rank_seq[0], rank_pos[0], 0     # before the reference start
            elif has_key[x]:
                # left flank: the largest reference rank strictly below the key
                # (a key equal to a rank lies parallel to that reference node)
                r0 = min(math.ceil(key[x]) - 1, last_rank)
                seq, rpos, on_ref = rank_seq[r0], rank_pos[r0] + length[rank_node[r0]], 0
            else:
                seq, rpos, on_ref = '-', '-', 0
            o.write(f'{names[x]}\t{new_id[x]}\t{flip[x]}\t{length[x]}\t{seg}\t{seq}\t{rpos}\t{on_ref}\n')

    bp_total = sum(length)
    removed_nodes = n - N
    bp_removed = sum(length[x] for x in range(n) if not keep[x])
    m_out = len(kept_edges)

    def pct(a, b):
        return round(100.0 * a / b, 4) if b else 0.0

    stats = {
        'input': {'nodes': n, 'links': m_in, 'bp': bp_total,
                  'haplotype_walks': n_walks, 'reference_paths': n_ref_paths},
        'orientation': {'nodes_flipped': sum(flip), 'pct_nodes_flipped': pct(sum(flip), n),
                        'walks_reversed': n_walks_reversed,
                        'walks_without_reference_node': n_walks_without_ref},
        'links_removed': {
            'inversion_links': n_inv, 'pct_inversion': pct(n_inv, m_in),
            'self_loops': n_self, 'pct_self_loops': pct(n_self, m_in),
            'merged_duplicates_after_orientation': n_dup,
            'cycle_breaking_back_edges': n_back, 'pct_back_edges': pct(n_back, m_in),
            'back_edges_on_reference': n_back_ref,
            'dropped_with_small_components': sum(1 for k in range(m_dir) if removed[k] == 2),
            'total_removed': m_in - m_out, 'pct_total_removed': pct(m_in - m_out, m_in),
        },
        'nodes_removed': {
            'components_before_filter': n_comp,
            'nodes': removed_nodes, 'pct_nodes': pct(removed_nodes, n),
            'bp': bp_removed, 'pct_bp': pct(bp_removed, bp_total),
        },
        'walks': {'transitions': trans, 'broken_transitions': broken,
                  'pct_broken': pct(broken, trans),
                  'haplotypes_with_a_broken_transition': len(haps_broken),
                  'batches': len(batch_order), 'unwalked_edges': len(unwalked)},
        'output': {'nodes': N, 'edges': m_out, 'sources': sources, 'sinks': sinks,
                   'acyclic': acyclic},
    }
    if args.pipeline_dir:
        write_pipeline_files(args, stats, batch_order, batches, ref_haps, unwalked,
                             eu, ev, new_id, length, time.time() - t_start)
        stats['pipeline_files'] = args.pipeline_dir

    with open(f'{args.out}.stats.json', 'w') as o:
        json.dump(stats, o, indent=2)

    lr = stats['links_removed']
    nr = stats['nodes_removed']
    w = stats['walks']
    print(f"""
Input:   {n:,} nodes, {m_in:,} links, {bp_total:,} bp, {n_walks} haplotype walks
Walks running against the reference (reverse-complemented contigs): {n_walks_reversed}; without any reference node: {n_walks_without_ref}
Links removed: {m_in - m_out:,} ({lr['pct_total_removed']}%)
  inversion links     {n_inv:,} ({lr['pct_inversion']}%)
  self-loops          {n_self:,} ({lr['pct_self_loops']}%)
  cycle back edges    {n_back:,} ({lr['pct_back_edges']}%), of which on reference: {n_back_ref}
  duplicates merged   {n_dup:,} (same directed edge after orientation; not lost information)
Nodes removed: {removed_nodes:,} ({nr['pct_nodes']}%), {bp_removed:,} bp ({nr['pct_bp']}%)
Walk transitions broken: {broken:,} of {trans:,} ({w['pct_broken']}%)
Output:  {N:,} nodes, {m_out:,} edges, {sources} source(s), {sinks} sink(s), acyclic={acyclic}
""")
    if not acyclic:
        sys.exit('error: output is not acyclic (this should not happen)')


if __name__ == '__main__':
    main()
