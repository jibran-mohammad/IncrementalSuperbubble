"""
superbubbleHistory.py

Turns the superbubble history written by the batch algorithm into tables a
reader can use: one row per superbubble lifespan, a summary per genome
(optionally with populations), and a BED file of the final superbubbles.

The history itself comes from the C++ program (perBatch mode, third argument
"history"; see lib/findSuperbubble/superbubbleHistory.h):

    ./test/incBSTBatch <dir>/incBSTBatch_time.txt perBatch history < numeric_edgelist.txt

=====================================================================
HOW TO RUN
=====================================================================

Run from the repository root.

    # minimal: lifespans and summary in input (edge list) vertex ids
    python3 Preprocessing/superbubbleHistory.py \
        --events  <dir>/incBSTBatch_history_events.tsv \
        --summary <dir>/incBSTBatch_history_summary.tsv \
        --out     <dir>/history

    # HPRC graph prepared by mcToDag.py: original node names, orientation,
    # reference coordinates and a BED file; with populations
    python3 Preprocessing/superbubbleHistory.py \
        --events  output/HPRC/chromosome_22/timeFile/perbatch/incBSTBatch_history_events.tsv \
        --summary output/HPRC/chromosome_22/timeFile/perbatch/incBSTBatch_history_summary.tsv \
        --nodemap input/HPRC/chromosome_22/Graphs/chromosome_22.nodemap.tsv \
        --metadata samples.tsv \
        --out     output/HPRC/chromosome_22/history/chr22

Options:
    --events     <ID>_history_events.tsv   (required)
    --summary    <ID>_history_summary.tsv  (required)
    --nodemap    <name>.nodemap.tsv of mcToDag.py (optional)
    --metadata   sample table, tab separated, with a header (optional):
                 a column "sample" and any other columns, e.g. "population",
                 "superpopulation"; a genome is matched to a sample if its
                 name starts with the sample id followed by '#', '.' or '_',
                 or contains it between '_' (e.g. GCA_..._NA19468_hap1_...)
    --out        output prefix (required)

=====================================================================
OUTPUT
=====================================================================

<out>_lifespans.tsv, one row per lifespan (a superbubble that dies and
appears again later has two rows):
    sb_id            running number
    entrance, exit   input vertex ids (edge list)
    entrance_node,   with --nodemap: original segment name with the
    exit_node        orientation in walk notation (>name forward, <name
                     reversed relative to the original graph)
    ref_seq, ref_start, ref_end
                     with --nodemap: reference coordinates (0-based, end
                     exclusive): from the entrance's position to the exit's
                     end (exit on the reference: position + length)
    n_inner          inner vertices when the superbubble appeared
    born_step, born_by, born_as
                     step (0 = base graph), genome, and how it appeared:
                     base / born (contains a vertex the genome added) /
                     promoted (only vertices that existed before)
    parent_at_birth  innermost superbubble containing it when it appeared
    died_step, died_by, absorbed_by
                     NA if it survives to the end; absorbed_by: the
                     innermost superbubble containing its former entrance
                     and exit after it died ("-": none)

<out>_summary.tsv: the summary of the program per genome plus
    sample and the --metadata columns (if given), cum_born, cum_promoted,
    cum_died (cumulative over the genomes, without the base graph)

<out>_final.bed (with --nodemap): the superbubbles alive at the end, BED9:
    chrom, start, end, name (entrance-exit), score (born_step, max 1000),
    strand ".", thickStart, thickEnd, itemRgb: colour by born_step (blue =
    base graph / early, red = late)

The script also prints a check: the superbubbles alive at the end (from the
lifespans) must equal the "alive" count of the last summary row.
"""
import argparse
import collections
import csv
import os
import sys


def read_nodemap(path):
    """new id -> (orig name, flipped, length, ref_seq, ref_pos, on_ref)"""
    nodes = {}
    with open(path) as f:
        header = f.readline().rstrip('\n').split('\t')
        col = {name: i for i, name in enumerate(header)}
        has_ref = 'ref_pos' in col
        for line in f:
            p = line.rstrip('\n').split('\t')
            new = int(p[col['new']])
            if new < 0:
                continue
            ref_seq = p[col['ref_seq']] if has_ref else '-'
            ref_pos = p[col['ref_pos']] if has_ref else '-'
            nodes[new] = (p[col['orig']], p[col['flipped']] == '1', int(p[col['length']]),
                          ref_seq, None if ref_pos == '-' else int(ref_pos),
                          has_ref and p[col['on_ref']] == '1')
    if not has_ref:
        print(f"[!] {path} has no reference coordinates (older mcToDag.py); "
              "run the preprocessing again for coordinates and the BED file", file=sys.stderr)
    return nodes, has_ref


def read_metadata(path):
    with open(path) as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    if not rows or 'sample' not in rows[0]:
        sys.exit(f"error: {path} needs a header with a column 'sample'")
    columns = [c for c in rows[0].keys() if c != 'sample']
    return {r['sample']: r for r in rows}, columns


def match_sample(genome, samples):
    """sample id of a genome name, or None"""
    for sep in ('#', '.', '_'):
        head = genome.split(sep)[0]
        if head in samples:
            return head
    for part in genome.split('_'):
        if part in samples:
            return part
    return None


def main():
    ap = argparse.ArgumentParser(description="Tables from a superbubble history")
    ap.add_argument('--events', required=True)
    ap.add_argument('--summary', required=True)
    ap.add_argument('--nodemap')
    ap.add_argument('--metadata')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    if os.path.dirname(args.out):
        os.makedirs(os.path.dirname(args.out), exist_ok=True)

    nodes, has_ref = read_nodemap(args.nodemap) if args.nodemap else ({}, False)

    def node_name(vid):
        n = nodes.get(int(vid) - 2)          # edge list id = new id + 2
        if n is None:
            return '-'
        return ('<' if n[1] else '>') + n[0]

    def coords(entrance, exit_):
        a, b = nodes.get(int(entrance) - 2), nodes.get(int(exit_) - 2)
        if a is None or b is None or a[4] is None or b[4] is None or a[3] != b[3]:
            return '-', '-', '-'
        start = a[4]
        end = b[4] + (b[2] if b[5] else 0)
        if end < start:
            start, end = end, start
        return a[3], start, max(end, start + 1)

    # ---- lifespans from the event log
    lifespans = []                     # dicts in order of appearance
    open_span = {}                     # (entrance, exit) -> index of the open lifespan
    last_step = 0
    with open(args.events) as f:
        next(f)
        for line in f:
            step, genome, ev, s, t, inner, parent = line.rstrip('\n').split('\t')
            step = int(step)
            last_step = max(last_step, step)
            k = (s, t)
            if ev == 'died':
                i = open_span.pop(k, None)
                if i is None:
                    sys.exit(f"error: {s}-{t} dies at step {step} but is not alive (inconsistent event log)")
                lifespans[i].update(died_step=step, died_by=genome, absorbed_by=parent)
            else:
                if k in open_span:
                    sys.exit(f"error: {s}-{t} appears at step {step} while alive (inconsistent event log)")
                open_span[k] = len(lifespans)
                lifespans.append(dict(entrance=s, exit=t, n_inner=inner, born_step=step, born_by=genome,
                                      born_as=ev, parent_at_birth=parent, died_step='NA', died_by='NA',
                                      absorbed_by='NA'))

    with open(args.out + '_lifespans.tsv', 'w') as o:
        o.write('sb_id\tentrance\texit\tentrance_node\texit_node\tref_seq\tref_start\tref_end\tn_inner\t'
                'born_step\tborn_by\tborn_as\tparent_at_birth\tdied_step\tdied_by\tabsorbed_by\n')
        for i, L in enumerate(lifespans, 1):
            seq, start, end = coords(L['entrance'], L['exit']) if nodes else ('-', '-', '-')
            o.write(f"{i}\t{L['entrance']}\t{L['exit']}\t{node_name(L['entrance']) if nodes else '-'}\t"
                    f"{node_name(L['exit']) if nodes else '-'}\t{seq}\t{start}\t{end}\t{L['n_inner']}\t"
                    f"{L['born_step']}\t{L['born_by']}\t{L['born_as']}\t{L['parent_at_birth']}\t"
                    f"{L['died_step']}\t{L['died_by']}\t{L['absorbed_by']}\n")

    # ---- summary with populations and cumulative counts
    meta, meta_cols = read_metadata(args.metadata) if args.metadata else ({}, [])
    unmatched = []
    with open(args.summary) as f, open(args.out + '_summary.tsv', 'w') as o:
        reader = csv.DictReader(f, delimiter='\t')
        cols = reader.fieldnames + (['sample'] + meta_cols if args.metadata else []) + \
            ['cum_born', 'cum_promoted', 'cum_died']
        o.write('\t'.join(cols) + '\n')
        cb = cp = cd = 0
        last_alive = None
        for r in reader:
            if int(r['step']) > 0:
                cb += int(r['born']); cp += int(r['promoted']); cd += int(r['died'])
            row = [r[c] for c in reader.fieldnames]
            if args.metadata:
                sample = match_sample(r['genome'], meta) if int(r['step']) > 0 else None
                if sample is None and int(r['step']) > 0:
                    unmatched.append(r['genome'])
                row += [sample or '-'] + [meta[sample][c] if sample else '-' for c in meta_cols]
            row += [str(cb), str(cp), str(cd)]
            o.write('\t'.join(row) + '\n')
            last_alive = int(r['alive'])

    # ---- BED of the final superbubbles
    alive_end = [L for L in lifespans if L['died_step'] == 'NA']
    if nodes and has_ref:
        def colour(step):
            f = step / last_step if last_step else 0.0
            return f"{int(255 * f)},0,{int(255 * (1 - f))}"
        n_bed = 0
        with open(args.out + '_final.bed', 'w') as o:
            o.write(f'track name="superbubbles" description="superbubbles of the final graph, '
                    f'colour = step of appearance (blue early, red late)" itemRgb="On"\n')
            rows = []
            for L in alive_end:
                seq, start, end = coords(L['entrance'], L['exit'])
                if seq == '-':
                    continue
                rows.append((seq, start, end, f"{L['entrance']}-{L['exit']}", L['born_step']))
            for seq, start, end, name, born in sorted(rows):
                o.write(f"{seq}\t{start}\t{end}\t{name}\t{min(born, 1000)}\t.\t{start}\t{end}\t{colour(born)}\n")
                n_bed += 1
        print(f"BED: {n_bed} superbubbles with coordinates -> {args.out}_final.bed")

    # ---- check and short report
    born_as = collections.Counter(L['born_as'] for L in lifespans)
    print(f"lifespans: {len(lifespans)} ({born_as['base']} base, {born_as['born']} born, "
          f"{born_as['promoted']} promoted), alive at the end: {len(alive_end)}")
    ok = (last_alive == len(alive_end))
    print(f"check: alive at the end = last summary row ({last_alive}): {'OK' if ok else 'MISMATCH'}")
    if args.metadata:
        print(f"genomes without a sample in {args.metadata}: {len(unmatched)}"
              + (f" (e.g. {', '.join(unmatched[:3])})" if unmatched else ""))
    print(f"written: {args.out}_lifespans.tsv, {args.out}_summary.tsv")
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
