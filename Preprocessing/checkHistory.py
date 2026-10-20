"""
checkHistory.py

Checks a superbubble history written by the batch algorithm
(incBSTBatch ... perBatch history; see lib/findSuperbubble/superbubbleHistory.h):
the event log is replayed from the base graph on, and
  - after every step the number of alive superbubbles must equal the count in
    the time file (SB_Count, which includes the artificial <0 1>),
  - no superbubble may die while not alive or be born while alive,
  - after the last step the alive set must equal the final superbubble list
    (<ID>_superbubbles.txt, which also contains <0 1>).

Usage (from the repository root):
    python3 Preprocessing/checkHistory.py \
        --events output/Human/chromosome_21/timeFile/perbatch/incBSTBatch_history_events.tsv \
        --time   output/Human/chromosome_21/timeFile/perbatch/incBSTBatch_time.txt \
        --final  output/Human/chromosome_21/timeFile/perbatch/incBSTBatch_superbubbles.txt

Exit code 0 if everything matches, 1 otherwise.
"""
import argparse
import collections
import sys


def main():
    ap = argparse.ArgumentParser(description="Replay a superbubble history and check it")
    ap.add_argument("--events", required=True, help="<ID>_history_events.tsv")
    ap.add_argument("--time", required=True, help="<ID>_time.txt of the same run")
    ap.add_argument("--final", required=True, help="<ID>_superbubbles.txt of the same run")
    args = ap.parse_args()

    per_step = collections.defaultdict(list)
    with open(args.events) as f:
        next(f)
        for line in f:
            step, genome, ev, s, t, inner, parent = line.rstrip("\n").split("\t")
            per_step[int(step)].append((ev, (s, t)))

    # superbubble count per step from the time file, without <0 1>
    counts = {}
    with open(args.time) as f:
        for line in f:
            p = line.split()
            if p and p[0] == "BaseGraph_Build":
                counts[0] = int(p[2]) - 1
            elif p and p[0].isdigit() and len(p) >= 10:
                counts[int(p[0])] = int(p[9]) - 1

    alive = set()
    bad_counts, died_not_alive, born_alive = [], 0, 0
    last_step = max(list(per_step) + list(counts))
    for step in range(0, last_step + 1):
        for ev, k in per_step.get(step, []):
            if ev == "died":
                if k not in alive:
                    died_not_alive += 1
                alive.discard(k)
            else:
                if k in alive:
                    born_alive += 1
                alive.add(k)
        if step in counts and len(alive) != counts[step]:
            bad_counts.append((step, len(alive), counts[step]))

    final = set()
    with open(args.final) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            s, t = line.split()
            if (s, t) != ("0", "1"):
                final.add((s, t))

    ok = not bad_counts and died_not_alive == 0 and born_alive == 0 and alive == final
    print(f"steps checked: {last_step + 1}")
    print(f"alive count = time-file count after every step: {not bad_counts}"
          + (f" (first differences, step/replayed/time file: {bad_counts[:3]})" if bad_counts else ""))
    print(f"died while not alive: {died_not_alive}, born while alive: {born_alive}")
    print(f"final set: replayed {len(alive)}, list {len(final)}, only replayed {len(alive - final)}, "
          f"only in list {len(final - alive)}")
    print("OK" if ok else "MISMATCH")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
