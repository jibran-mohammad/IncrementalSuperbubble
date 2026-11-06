#!/usr/bin/env python3
"""
superbubbleAlleles.py - variant type, alleles and edge support of the
superbubbles of the final graph, from the haplotype walks of the original
Minigraph-Cactus GFA (HPRC). Complements the lifespans of
Preprocessing/superbubbleHistory.py (one row per final superbubble, joined
on entrance / exit).

=====================================================================
HOW TO RUN
=====================================================================

Run from the repository root, after superbubbleHistory.py (its lifespans
give the final superbubbles):

    python3 Preprocessing/superbubbleAlleles.py \\
        --gfa       input/HPRC/chromosome_22/Graphs/chromosome_22.gfa[.gz] \\
        --nodemap   input/HPRC/chromosome_22/Graphs/chromosome_22.nodemap.tsv \\
        --graph     input/HPRC/chromosome_22/Graphs/numeric_edgelist.txt \\
        --lifespans output/HPRC/chromosome_22/history/chr22_lifespans.tsv \\
        --ref GRCh38 \\
        --out       output/HPRC/chromosome_22/history/chr22

Options:
    --gfa        the original GFA (W or P lines with the haplotype walks,
                 S lines with the sequences), as given to mcToDag.py
    --nodemap    <name>.nodemap.tsv of mcToDag.py (original name -> new id,
                 flipped)
    --graph      numeric_edgelist.txt (the final DAG)
    --lifespans  <out>_lifespans.tsv of superbubbleHistory.py
    --ref        reference sample in the walks (default: GRCh38)
    --out        output prefix

=====================================================================
OUTPUT: <out>_alleles.tsv, one row per superbubble of the final graph
=====================================================================

    entrance, exit     edge-list ids (as in the lifespans)
    var_type           SNP, MNP, insertion, deletion, complex (every
                       non-reference allele against the reference allele:
                       ref 1 letter / other 1 letter = SNP; same length >= 2
                       = MNP; other empty = deletion; ref empty = insertion;
                       different lengths, both non-empty = complex); the
                       common type if all non-reference alleles agree, else
                       mixed; nested = contains other superbubbles (not
                       classified, as vg deconstruct's levels); nonref = the
                       reference does not cross the site; none = every
                       haplotype carries the reference allele
    n_alleles          distinct allele sequences, the reference included
                       (> 2: multi-allelic); an allele is the sequence
                       between entrance and exit (both excluded) spelled by a
                       walk that crosses the site completely
    n_hap_traversing   haplotypes (sample#hap) crossing the site completely
                       (a haplotype split into several contigs counts once)
    alt_count          haplotypes among them whose allele differs from the
                       reference allele (the reference itself not counted)
    min_edge_support   the smallest number of haplotypes using an edge inside
                       the site (entrance .. exit); 1 = a singleton edge

Derived: allele frequency = alt_count / (n_hap_traversing - 1, without the
reference); missing = haplotypes in the graph - n_hap_traversing.

How a walk is followed: its steps are mapped to the DAG (new ids, the
orientation of mcToDag.py); a walk that runs mostly against the DAG
direction (a reverse-complemented contig) is reversed. At the entrance of a
final superbubble the walk is followed while it stays strictly inside the
site (positions between entrance and exit in a DFS order of the final DAG,
in which every superbubble is that range); reaching the exit completes the
allele, leaving the range or a node removed by the DAG normalization
abandons it.
"""
import argparse
import collections
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mcToDag import opener, path_steps                          # noqa: E402
from superbubbleHistory import read_edges, dfs_order            # noqa: E402

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def classify(ref, alt):
    """type of one non-reference allele against the reference allele"""
    if len(ref) == 1 and len(alt) == 1:
        return "SNP"
    if len(alt) == 0:
        return "deletion"
    if len(ref) == 0:
        return "insertion"
    if len(ref) == len(alt):
        return "MNP"
    return "complex"


def main():
    ap = argparse.ArgumentParser(description="Variant type and alleles of the final superbubbles")
    ap.add_argument("--gfa", required=True)
    ap.add_argument("--nodemap", required=True)
    ap.add_argument("--graph", required=True)
    ap.add_argument("--lifespans", required=True)
    ap.add_argument("--ref", default="GRCh38")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    # ---- final DAG: DFS order (every superbubble = range of positions) and successors
    edges = read_edges(args.graph)
    pos, order = dfs_order(edges)
    succ = collections.defaultdict(list)
    for u, v in edges:
        succ[u].append(v)

    # ---- final superbubbles from the lifespans
    sb_of = {}                     # entrance id -> (exit id, a, b)
    ranges = []
    with open(args.lifespans) as f:
        head = f.readline().rstrip("\n").split("\t")
        c = {n: i for i, n in enumerate(head)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[c["died_step"]] != "NA":
                continue
            s, t = int(p[c["entrance"]]), int(p[c["exit"]])
            a, b = pos[s], pos[t]
            if a < 0 or b < 0 or a >= b:
                continue
            sb_of[s] = (t, a, b)
            ranges.append((a, b, s))
    # superbubbles containing others (nested): a range with another one inside
    ranges.sort()
    has_child = set()
    stack = []
    for a, b, s in ranges:
        while stack and stack[-1][0] <= a:
            stack.pop()
        if stack:
            has_child.add(stack[-1][1])
        stack.append((b, s))
    print(f"final superbubbles: {len(sb_of)} ({len(has_child)} with nested superbubbles)")

    # ---- node map: original name -> (edge-list id, flipped)
    node = {}
    with open(args.nodemap) as f:
        head = f.readline().rstrip("\n").split("\t")
        c = {n: i for i, n in enumerate(head)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            new = int(p[c["new"]])
            node[p[c["orig"]]] = (new + 2 if new >= 0 else -1, p[c["flipped"]] == "1")

    # ---- pass 1: sequences of the nodes inside final superbubbles (DAG orientation)
    inside = bytearray(len(pos))
    for s, (t, a, b) in sb_of.items():
        for q in range(a + 1, b):
            inside[order[q]] = 1
    seq = {}
    with opener(args.gfa) as fh:
        for line in fh:
            if line[0] != "S":
                continue
            f = line.rstrip("\n").split("\t", 3)
            x, flipped = node.get(f[1], (-1, False))
            if x >= 0 and x < len(inside) and inside[x]:
                s_ = f[2] if f[2] != "*" else ""
                seq[x] = s_.translate(COMP)[::-1] if flipped else s_
    print(f"sequences of {len(seq)} inner nodes read")

    # ---- pass 2: walks -> alleles per superbubble and haplotypes per edge
    alleles = collections.defaultdict(dict)        # entrance -> {haplotype: allele}
    ref_allele = {}
    edge_haps = collections.defaultdict(set)       # (u, v) inside a superbubble -> haplotypes
    hap_ids = {}
    n_walks = n_rev = 0
    with opener(args.gfa) as fh:
        for line in fh:
            if line[0] not in "WP":
                continue
            f = line.rstrip("\n").split("\t")
            sample, hap, steps = path_steps(f)
            n_walks += 1
            h = hap_ids.setdefault(hap, len(hap_ids))
            is_ref = sample == args.ref
            mapped = [(node.get(nm, (-1, False))[0], rv != node.get(nm, (-1, False))[1]) for nm, rv in steps]
            # mapped: (id, backward against the DAG); a walk mostly backward is reversed
            back = sum(1 for x, bw in mapped if x >= 0 and bw)
            if back * 2 > sum(1 for x, _ in mapped if x >= 0):
                mapped = [(x, not bw) for x, bw in reversed(mapped)]
                n_rev += 1
            w = [x for x, _ in mapped]
            for i in range(len(w) - 1):
                u, v = w[i], w[i + 1]
                if u >= 0 and v >= 0 and inside_or_end(u, v, pos, sb_of, inside):
                    edge_haps[(u, v)].add(h)
            for i, x in enumerate(w):
                if x < 0 or x not in sb_of:
                    continue
                t, a, b = sb_of[x]
                parts = []
                j = i + 1
                complete = False
                while j < len(w):
                    y = w[j]
                    if y == t:
                        complete = True
                        break
                    if y < 0 or not (a < pos[y] < b):
                        break
                    parts.append(seq.get(y, ""))
                    j += 1
                if not complete:
                    continue
                allele = "".join(parts)
                if is_ref:
                    ref_allele.setdefault(x, allele)
                else:
                    alleles[x].setdefault(h, allele)
    print(f"walks: {n_walks} ({n_rev} reversed), haplotypes: {len(hap_ids)}")

    # ---- edge support: smallest number of haplotypes on an edge inside the site
    def min_support(s, t, a, b):
        m = None
        for q in range(a, b):
            x = order[q]
            for y in succ.get(x, ()):
                if pos[y] <= b:
                    k = len(edge_haps.get((x, y), ()))
                    m = k if m is None else min(m, k)
        return m if m is not None else 0

    types = collections.Counter()
    with open(args.out + "_alleles.tsv", "w") as o:
        o.write("entrance\texit\tvar_type\tn_alleles\tn_hap_traversing\talt_count\tmin_edge_support\n")
        for s in sorted(sb_of):
            t, a, b = sb_of[s]
            haps = alleles.get(s, {})
            ref = ref_allele.get(s)
            n_trav = len(haps) + (1 if ref is not None else 0)
            distinct = set(haps.values()) | ({ref} if ref is not None else set())
            if s in has_child:
                vt = "nested"
            elif ref is None:
                vt = "nonref"
            else:
                kinds = {classify(ref, al) for al in set(haps.values()) if al != ref}
                vt = "none" if not kinds else (kinds.pop() if len(kinds) == 1 else "mixed")
            alt = sum(1 for al in haps.values() if ref is not None and al != ref)
            types[vt] += 1
            o.write(f"{s}\t{t}\t{vt}\t{len(distinct)}\t{n_trav}\t{alt if ref is not None else 'NA'}\t"
                    f"{min_support(s, t, a, b)}\n")
    total = sum(types.values())
    print("var_type: " + ", ".join(f"{k} {v} ({100.0 * v / total:.1f}%)" for k, v in types.most_common()))
    print(f"written: {args.out}_alleles.tsv")


def inside_or_end(u, v, pos, sb_of, inside):
    """an edge that can lie inside a final superbubble: one end inside one, or
    from an entrance (only these edges are needed for the edge support)"""
    return inside[u] or inside[v] or u in sb_of


if __name__ == "__main__":
    main()
