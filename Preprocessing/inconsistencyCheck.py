#!/usr/bin/env python3
"""
gfa_check_consistency.py — sanity-check a GFA for link/walk inconsistencies.

Checks:
  1) Missing links: oriented adjacencies present in W/P but absent in L.
  2) Orphan links: oriented adjacencies present in L but never traversed by W/P.
  3) Node coverage: S nodes not referenced by any W/P; W/P/L referencing unknown S.

Exit status:
  0 if no issues; 1 if any inconsistency found.

Usage:
  python3 gfa_check_consistency.py path/to/graph.gfa [--max-examples 50] [--report-prefix reports/out]

Outputs (optional with --report-prefix PREFIX):
  PREFIX_missing_in_L.tsv
  PREFIX_orphan_L.tsv
  PREFIX_nodes_not_in_walks.tsv
  PREFIX_refs_to_unknown_S.tsv
"""

import argparse
import re
import sys
from pathlib import Path

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("gfa", type=Path, help="Input GFA file")
    ap.add_argument("--max-examples", type=int, default=50, help="Max examples to print per category")
    ap.add_argument("--report-prefix", type=Path, default=None, help="If set, write TSV reports with this prefix")
    return ap.parse_args()

# Regex helpers
W_TOKEN_RE = re.compile(r'([<>])(s[^\s><,]+)')      # >s123 or <s456 (contiguous)
P_TOKEN_RE = re.compile(r'(s[^\s,]+?)([+-])')       # s123+ or s456-
TAB = "\t"

def main():
    args = parse_args()
    if not args.gfa.exists():
        print(f"ERROR: File not found: {args.gfa}", file=sys.stderr)
        sys.exit(2)

    # Collectors
    S_nodes = set()
    link_edges = set()    # tuples: (a, a_or, b, b_or)
    walk_edges = set()    # edges implied by W/P
    nodes_in_walks = set()
    refs_to_unknown_S = set()  # any segment referenced by W/P/L without an S line
    malformed_lines = 0

    # First pass: read all lines
    with args.gfa.open("r") as fh:
        for ln, line in enumerate(fh, 1):
            if not line or line[0] == '#':
                continue

            if line.startswith('S\t'):
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 2:
                    S_nodes.add(parts[1])
                else:
                    malformed_lines += 1

            elif line.startswith('L\t'):
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 5:
                    a, a_or, b, b_or = parts[1], parts[2], parts[3], parts[4]
                    link_edges.add((a, a_or, b, b_or))
                    # Track unknown references for later reporting
                    if a not in S_nodes:
                        refs_to_unknown_S.add(("L", a))
                    if b not in S_nodes:
                        refs_to_unknown_S.add(("L", b))
                else:
                    malformed_lines += 1

            elif line.startswith('W\t'):
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 7:
                    walk_field = parts[6]
                    toks = W_TOKEN_RE.findall(walk_field)  # [('>', 's1'), ('<','s2'), ...]
                    # node presence
                    for _, seg in toks:
                        nodes_in_walks.add(seg)
                        if seg not in S_nodes:
                            refs_to_unknown_S.add(("W", seg))
                    # oriented adjacencies
                    for (o1, a), (o2, b) in zip(toks, toks[1:]):
                        ao = '+' if o1 == '>' else '-'
                        bo = '+' if o2 == '>' else '-'
                        walk_edges.add((a, ao, b, bo))
                else:
                    malformed_lines += 1

            elif line.startswith('P\t'):
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 3:
                    path_field = parts[2]  # s123+,s124-,...
                    toks = P_TOKEN_RE.findall(path_field)  # [('s123','+'), ...]
                    for seg, sign in toks:
                        nodes_in_walks.add(seg)
                        if seg not in S_nodes:
                            refs_to_unknown_S.add(("P", seg))
                    for (a, s1), (b, s2) in zip(toks, toks[1:]):
                        ao = '+' if s1 == '+' else '-'
                        bo = '+' if s2 == '+' else '-'
                        walk_edges.add((a, ao, b, bo))
                else:
                    malformed_lines += 1

            else:
                # other lines (H, a, etc.) ignored for this check
                pass

    # Compute inconsistencies
    missing_in_L = walk_edges - link_edges
    orphan_L = link_edges - walk_edges
    nodes_not_in_walks = S_nodes - nodes_in_walks

    # Summary
    print("=== GFA Consistency Check ===")
    print(f"S nodes                      : {len(S_nodes)}")
    print(f"W/P implied oriented edges   : {len(walk_edges)}")
    print(f"L oriented edges             : {len(link_edges)}")
    print(f"Missing links (W∖L)          : {len(missing_in_L)}")
    print(f"Orphan links (L∖W)           : {len(orphan_L)}")
    print(f"S nodes not in any walk      : {len(nodes_not_in_walks)}")
    print(f"Refs to unknown S (W/P/L)    : {len(refs_to_unknown_S)}")
    if malformed_lines:
        print(f"Malformed lines (skipped)    : {malformed_lines}")

    maxn = args.max_examples

    if missing_in_L:
        print("\n-- Examples: Missing oriented links (present in W/P, absent in L) --")
        for i, (a, ao, b, bo) in enumerate(sorted(missing_in_L)):
            if i >= maxn: break
            print(f"{a}{ao} -> {b}{bo}")

    if orphan_L:
        print("\n-- Examples: Orphan oriented links (present in L, never in W/P) --")
        for i, (a, ao, b, bo) in enumerate(sorted(orphan_L)):
            if i >= maxn: break
            print(f"{a}{ao} -> {b}{bo}")

    if nodes_not_in_walks:
        print("\n-- Examples: S nodes not seen in any walk --")
        for i, s in enumerate(sorted(nodes_not_in_walks)):
            if i >= maxn: break
            print(s)

    if refs_to_unknown_S:
        print("\n-- Examples: Walk/Link references to missing S nodes --")
        for i, (src, seg) in enumerate(sorted(refs_to_unknown_S)):
            if i >= maxn: break
            print(f"{src}: {seg}")

    # Optional TSV reports
    if args.report_prefix:
        prefix = args.report_prefix
        prefix.parent.mkdir(parents=True, exist_ok=True)

        def write_tsv(path, rows, header):
            with path.open("w") as fh:
                if header: fh.write(header + "\n")
                for r in rows:
                    if isinstance(r, (tuple, list)):
                        fh.write(TAB.join(map(str, r)) + "\n")
                    else:
                        fh.write(str(r) + "\n")

        write_tsv(prefix.with_name(prefix.name + "_missing_in_L.tsv"),
                  [(a, ao, b, bo) for (a, ao, b, bo) in sorted(missing_in_L)],
                  "segA\toriA\tsegB\toriB")

        write_tsv(prefix.with_name(prefix.name + "_orphan_L.tsv"),
                  [(a, ao, b, bo) for (a, ao, b, bo) in sorted(orphan_L)],
                  "segA\toriA\tsegB\toriB")

        write_tsv(prefix.with_name(prefix.name + "_nodes_not_in_walks.tsv"),
                  sorted(nodes_not_in_walks),
                  "segment")

        write_tsv(prefix.with_name(prefix.name + "_refs_to_unknown_S.tsv"),
                  [(src, seg) for (src, seg) in sorted(refs_to_unknown_S)],
                  "source\tsegment")

        print(f"\nReports written with prefix: {prefix}")

    # Decide exit code
    has_issues = bool(missing_in_L or orphan_L or nodes_not_in_walks or refs_to_unknown_S)
    sys.exit(1 if has_issues else 0)

if __name__ == "__main__":
    main()
