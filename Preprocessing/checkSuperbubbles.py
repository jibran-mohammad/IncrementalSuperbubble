r"""
checkSuperbubbles.py

Checks the superbubbles that the C++ algorithms report against an
independent reference, after every update: the reference recomputes all
superbubbles of the current graph from scratch with the definition-based
algorithm of Onodera et al. (2013), on the graph the algorithms work on
(base edges + the updates so far, plus the artificial source 0 with an
edge to every vertex without a parent, the artificial sink 1 with an edge
from every vertex without a child, and the edge 0 -> 1).

=====================================================================
HOW TO RUN
=====================================================================

Run from the repository root, with the virtual environment active
(source venv/bin/activate), after the algorithm has written its timing
file (scripts/superbubbleRunner.py, or a binary run by hand).

    # E. coli, per-edge results of all algorithms
    python3 Preprocessing/checkSuperbubbles.py \
        --edgelist input/Ecoli/Graphs/numeric_edgelist.txt -m perEdge \
        --time output/Ecoli/timeFile/perEdge/static_time.txt \
               output/Ecoli/timeFile/perEdge/incMinIBST_time.txt

    # chromosome 21, per-haplotype results
    python3 Preprocessing/checkSuperbubbles.py \
        --edgelist input/Human/chromosome_21/Graphs/numeric_edgelist.txt -m perBatch \
        --time output/Human/chromosome_21/timeFile/perbatch/static_time.txt

    # a binary run by hand: ./test/static my_time.txt perEdge < <edge list>
    python3 Preprocessing/checkSuperbubbles.py --edgelist <edge list> -m perEdge --time my_time.txt

Options:

    --edgelist   numeric_edgelist.txt the algorithm was run on      (required)
    -m           perEdge | perBatch, as the algorithm was run         (required)
    --time       one or more <ID>_time.txt timing files               (required)
    --every N    check only every N-th update (the reference is slow
                 on large graphs; the last update is always checked)  (default: 1)

=====================================================================
OUTPUT
=====================================================================

For every timing file: the number of checked updates whose superbubble
count differs from the reference, the first such update, and, if
<ID>_superbubbles.txt is next to the timing file, a comparison of the
final graph's superbubbles (pairs only in the algorithm / only in the
reference). Exit code 1 if anything differs.

Counts: the algorithms report the pair <0 1> (artificial source and
sink) and no other pair with an artificial vertex, so their count minus 1
is compared with the reference's number of superbubbles without 0 or 1.
"""

import argparse
import collections
import pathlib
import sys


def read_edge_list(path):
    """Base edges and, per haplotype, its edges (numeric_edgelist.txt)."""

    lines = [l.split() for l in open(path)]

    # header "<base edges> <nodes>": the first line with two numbers
    i = next(k for k, x in enumerate(lines) if len(x) == 2 and x[0].isdigit() and x[1].isdigit())
    m = int(lines[i][0])
    base = [tuple(map(int, x)) for x in lines[i + 1:i + 1 + m]]

    haplotypes = []
    k = i + 1 + m
    while k < len(lines):
        x = lines[k]
        # "# HAPLOTYPE: <name>" (separator lines "# HAPLOTYPE: ####" are skipped)
        if len(x) >= 3 and x[0] == '#' and x[1] == 'HAPLOTYPE:' and not x[2].startswith('#'):
            k += 1
            while k < len(lines) and not lines[k]:
                k += 1
            # a genome that adds no new edge has no count line (and no edges)
            if k >= len(lines) or not lines[k][0].isdigit():
                haplotypes.append([])
                continue
            count = int(lines[k][0])
            haplotypes.append([tuple(map(int, e)) for e in lines[k + 1:k + 1 + count]])
            k += 1 + count
        else:
            k += 1

    return base, haplotypes


def reference_superbubbles(out, inn):
    """
    All superbubbles (s, t) of the graph given by out/inn (artificial source
    and sink already added), with the algorithm of Onodera et al. (2013):
    from every entrance candidate s, traverse until a unique exit is found.
    """

    result = set()

    for s in list(out):
        if not out[s]:
            continue
        seen = {s}
        visited = set()
        stack = [s]
        while stack:
            v = stack.pop()
            visited.add(v)
            seen.discard(v)
            if not out[v]:
                break                       # tip: no superbubble from s
            cycle = False
            for u in out[v]:
                if u == s:
                    cycle = True            # cycle through s
                    break
                seen.add(u)
                if inn[u] <= visited:       # all parents of u visited
                    stack.append(u)
            if cycle:
                break
            if len(stack) == 1 and seen == {stack[0]}:
                t = stack[0]
                if s not in out[t]:
                    result.add((s, t))
                break

    return result


class Graph:
    """Real edges plus the artificial source 0 / sink 1 as the algorithms define them."""

    def __init__(self):
        self.out = collections.defaultdict(set)
        self.inn = collections.defaultdict(set)

    def add(self, u, v):
        self.out[u].add(v)
        self.inn[v].add(u)
        self.out[v]
        self.inn[u]

    def augmented(self):
        out = {v: set(w) for v, w in self.out.items()}
        inn = {v: set(w) for v, w in self.inn.items()}
        out.setdefault(0, set()); inn.setdefault(0, set())
        out.setdefault(1, set()); inn.setdefault(1, set())
        for v in list(self.out):
            if not self.inn[v]:
                out[0].add(v); inn[v].add(0)
            if not self.out[v]:
                out[v].add(1); inn[1].add(v)
        out[0].add(1); inn[1].add(0)
        return out, inn


def real(pairs):
    """Superbubbles without the artificial vertices 0 and 1."""
    return {(s, t) for s, t in pairs if s > 1 and t > 1}


def read_counts(time_file, haplotypes, mode):
    """
    Superbubble count after each update (column SB_Count) of a timing file:
    one per edge (perEdge) or per genome (perBatch).

    The rows are read in step with the genomes of the edge list: a genome
    with edges has one row per edge (perEdge) or one row (perBatch); a
    genome that adds no edge has one short row "<k> <name> 0 <SB_Count>"
    in both modes, which is an update only in perBatch mode.
    """

    lines = open(time_file).read().splitlines()
    header = lines[2].split()
    # Haplotype names may contain spaces, so SB_Count is located from the
    # end of the row: the number of columns after it is the same in every row.
    after = len(header) - header.index("SB_Count") - 1
    # update rows start with a number and have at least 4 fields (the lines
    # "<edges> <vertices>" after each genome have 2)
    rows = iter(l.split() for l in lines[3:] if len(l.split()) >= 4 and l.split()[0].isdigit())

    counts = []
    for h in haplotypes:
        if not h:
            r = next(rows)
            if mode == "perBatch":
                counts.append(int(r[-1]))
        else:
            for _ in range(len(h) if mode == "perEdge" else 1):
                counts.append(int(next(rows)[-1 - after]))
    return counts


def read_final_list(time_file):
    path = pathlib.Path(str(time_file).replace("_time.txt", "_superbubbles.txt"))
    if path == pathlib.Path(time_file) or not path.exists():
        return None
    return {tuple(map(int, l.split()[:2])) for l in open(path) if l.strip() and not l.startswith('#')}


def main():

    parser = argparse.ArgumentParser(description="Check superbubble counts of the algorithms "
                                     "after every update against an independent reference.")
    parser.add_argument("--edgelist", required=True)
    parser.add_argument("-m", "--mode", choices=["perEdge", "perBatch"], required=True)
    parser.add_argument("--time", nargs="+", required=True)
    parser.add_argument("--every", type=int, default=1)
    args = parser.parse_args()

    base, haplotypes = read_edge_list(args.edgelist)

    # the updates: single edges (perEdge) or a haplotype's edges (perBatch)
    if args.mode == "perEdge":
        updates = [[e] for h in haplotypes for e in h]
    else:
        updates = [h for h in haplotypes]

    algorithms = {t: read_counts(t, haplotypes, args.mode) for t in args.time}
    for t, c in algorithms.items():
        if len(c) != len(updates):
            print(f"[!] {t}: {len(c)} updates in the timing file, {len(updates)} in the edge list "
                  f"(wrong mode or edge list?)")

    g = Graph()
    for u, v in base:
        g.add(u, v)

    differ = {t: [] for t in args.time}
    checked = 0
    final = None

    for k, edges in enumerate(updates):
        for u, v in edges:
            g.add(u, v)
        last = (k == len(updates) - 1)
        if k % args.every and not last:
            continue
        ref = real(reference_superbubbles(*g.augmented()))
        checked += 1
        if last:
            final = ref
        for t, counts in algorithms.items():
            if k < len(counts) and counts[k] - 1 != len(ref):
                differ[t].append((k + 1, counts[k] - 1, len(ref)))
        if checked % 200 == 0:
            print(f"  ... {k + 1}/{len(updates)} updates", file=sys.stderr)

    failed = False
    print(f"{args.edgelist} ({args.mode}): {len(updates)} updates, {checked} checked against the reference")
    for t in args.time:
        d = differ[t]
        if d:
            failed = True
            first = d[0]
            print(f"  {t}: count differs on {len(d)} of {checked} checked updates; "
                  f"first: update {first[0]} (algorithm {first[1]}, reference {first[2]})")
        else:
            print(f"  {t}: count matches the reference on all {checked} checked updates")
        listed = read_final_list(t)
        if listed is not None and final is not None:
            mine = real(listed)
            only_a, only_r = mine - final, final - mine
            if only_a or only_r:
                failed = True
            print(f"     final graph: {len(mine & final)} in both, {len(only_a)} only in the algorithm"
                  f"{' ' + str(sorted(only_a)[:5]) if only_a else ''}, {len(only_r)} only in the reference"
                  f"{' ' + str(sorted(only_r)[:5]) if only_r else ''}")

    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
