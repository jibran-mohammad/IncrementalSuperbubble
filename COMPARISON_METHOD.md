# How the algorithm is compared with an external tool

This document describes, from the code, how the best algorithm of this
repository (`BEST_ALGORITHM` in `scripts/pipelineRunner.py`, currently
`incMinIBST` = Inc-MinI-BST) is compared with an external superbubble tool
(currently BubbleGun): what is timed on each side and why, how memory and
the superbubbles are compared, how to read the speedup, and which earlier
results are not comparable with current ones.

How the graphs and updates themselves are prepared is described in
[GRAPH_PREPARATION.md](GRAPH_PREPARATION.md); algorithm names in
[ALGORITHMS.md](ALGORITHMS.md). The commands are in the README
("Rerunning All Experiments" and Step 3c).


## 1. What is compared

Both sides process **the same sequence of graphs**: the base graph
(= the reference walk), then one update at a time until the final graph
(the union of all genome walks). The sequence is produced once by the
preprocessing (`pipelineRunner.py` / `graphRunner.py`):

- the algorithm reads `numeric_edgelist.txt`;
- the tool gets, after every update, a GFA with the base graph
  (`reduced_graph.gfa`), the edges added so far (`unique_elements.txt` in
  the order of `ordered_haplotypes.txt`) and an `S` line for every segment
  they introduce (`Preprocessing/graphUpdater.py`).

Two modes:

| Mode | One update | Number of updates |
|---|---|---|
| perEdge | one edge | number of added edges (E. coli 2,236; chromosome 21 6,863) |
| perBatch | all edges of one genome | number of genomes (E. coli 5; chromosome 21 63) |

A static tool such as BubbleGun recomputes all superbubbles of the whole
graph after every update; the dynamic algorithm keeps its data structures
and only processes the change.


## 2. Time

### Algorithm: `AlgorithmTime` (`scripts/pipelineRunner.py`)

- **Build:** `BENCHMARK=1` (`-O3`; no per-step timers, no debug output, no
  output files). Compilation is not timed.
- **What is timed:** the wall-clock time of the program, from start to
  exit: reading the edge list, building the base graph and computing its
  superbubbles, and processing all updates (in perEdge mode the
  superbubbles are recomputed after every edge, in perBatch after every
  genome). Program start-up is about 1 ms.
- **Repetitions:** `--repeat` runs (default 10); `AlgorithmTime` is the
  mean, `AlgorithmTimeStd` the standard deviation (E. coli perEdge:
  0.29 s +- 0.01 s; chromosome 21 perEdge: 4.75 s +- 0.13 s).

### Tool: `TotalTime` (`scripts/superbubbleComparator.py`, `Preprocessing/superbubbleTools/bubblegunTool.py`)

- **How it runs:** BubbleGun is used as a Python library inside the
  comparator, so Python start-up and importing BubbleGun (about 50 ms per
  run) are **not** timed: a user of a static tool on a changing graph
  would keep it loaded.
- **What is timed, per update:** reading the whole GFA into BubbleGun's
  graph + `find_bubbles` (finding all bubbles from scratch). This is what
  a static tool has to redo after every update. Building chains and
  nesting (`ChainsTime`) is extra BubbleGun work, recorded but **not**
  included. Writing the temporary GFA is done by the comparator and is not
  timed.
- `TotalTime` = base graph + all updates (as for the algorithm, which also
  includes its base graph). One run, no repetition; it is the sum of
  thousands of runs in perEdge mode.

### Tool: BubbleFinder (`Preprocessing/superbubbleTools/bubblefinderTool.py`)

- **How it runs:** BubbleFinder (C++/Rust, linear-time superbubbles via
  SPQR trees; Sena, Politov et al., arXiv:2511.21919; release 1.0.5) is a
  separate program, run once per update on the **same temporary GFA** as
  BubbleGun: `BubbleFinder superbubbles --directed --gfa-directed -j 1`.
  One thread, like our algorithms and BubbleGun.
- **What is timed, per update:** BubbleFinder's own timers
  (`--report-json`): reading the GFA (`io/read_graph`) + all phases of
  finding the superbubbles (`sb/...`). The same work as for BubbleGun;
  program start-up (about 10 ms) and writing the output file are **not**
  included. `TotalTime` as for BubbleGun.
- **Checked before use:** on the worked example (6 graphs: nesting of
  depth 2, separate components, superbubbles disappearing and
  reappearing) BubbleFinder finds exactly the reference superbubbles, and
  on chromosome 21's final graph the same 9,998 superbubbles as
  Inc-MinI-BST (2,573 nontrivial, as BubbleGun).

### Not included on either side

**Preprocessing** (GFA -> edge list / reduced GFA, `PreprocessingTime`)
produces the input of both, so it is reported once, as a separate column,
and added to neither time.

### Speedup

`Speedup = TotalTime (tool) / AlgorithmTime (algorithm)`, per dataset and
mode (`performance_<mode>.tex` of `Preprocessing/toolComparisonToLatex.py`).


## 3. Memory

| | Algorithm | Tool |
|---|---|---|
| Value | `PeakMemory(KB)` in `<ID>_time.txt` (`superbubbleRunner.py`, normal build, mean of 100 runs) | `PeakMemory(KB)` in `<tool>_time.txt` |
| Measured | peak resident memory of the program over the whole run (`VmHWM` in `/proc/self/status`, read just before exit) | BubbleGun: peak resident memory of a fresh Python process that reads the graph and finds the bubbles; BubbleFinder: peak resident memory of its process (from its `--report-json`); in a separate run at the base graph and after each genome, the largest value |
| Includes | the C++ runtime | BubbleGun: the Python interpreter + BubbleGun (`BaselineMemory(KB)`, about 14 MB, measured before the graph is read); BubbleFinder: the program and its libraries (`BaselineMemory(KB)` = its peak on a graph with one edge) |

Both are whole-process peaks, measured the same way (`VmHWM`). The
tool's memory is measured in separate, untimed runs, so measuring it does
not change the tool's time. On small graphs the tool's value is dominated
by the Python baseline; the graph-dependent part is `PeakMemory -
BaselineMemory`.


## 4. Superbubbles (agreement)

`agreement.tex`: the superbubbles of the **final graph** found by the
algorithm (`<ID>_superbubbles.txt`) and by the tool (last block of
`<tool>_superbubbles.txt`), compared as sets:

- pairs are written as GFA segment names (input id k = segment `s{k-1}`)
  and compared **unordered** (tools may report the ends in either order);
- left out on **both** sides, using the final graph from
  `numeric_edgelist.txt`: pairs with the artificial source/sink (ids 0 and
  1, which only the algorithm has) and **single-edge superbubbles** (the
  entrance's only out-neighbour is the exit), which BubbleGun does not
  report;
- BubbleGun's three bubble types (simple, super, insertion) all count.

If the algorithm is correct and the tool reports every superbubble, the
columns "Only <algorithm>" and "Only <tool>" are 0;
`agreement_details_<dataset>.txt` lists any differing superbubbles.

The agreement table checks only the final graph. Correctness after every
update is checked with `Preprocessing/checkSuperbubbles.py`, which
recomputes the superbubbles after every update with an independent
reference (definition-based algorithm of Onodera et al.) and compares the
count of each algorithm (column `SB_Count` of `<ID>_time.txt`): an error
can cancel out by the final graph (on 10.10.2026 the final graphs agreed,
but the static algorithm was wrong on 971 of 2,236 E. coli updates).


## 5. How to read the speedup

- **perEdge** is the case dynamic algorithms are made for: many small
  updates. The tool reads and processes the whole graph thousands of times;
  the algorithm processes one edge each time.
- **perBatch** has few, large updates (a whole genome each). The tool
  recomputes only once per genome while the algorithm still processes
  every edge of the genome, so the speedup is much smaller and can be
  about 1 on small datasets.
- Speedups of different modes or datasets are not comparable with each
  other; report them per mode, together with the number of updates and
  the time per update.

Reference values (10.10.2026, commit 17895c1; Inc-MinI-BST vs BubbleGun):

| Dataset, mode | Updates | BubbleGun per update | Inc-MinI-BST per update | Speedup (end-to-end) |
|---|---|---|---|---|
| E. coli, perEdge | 2,236 | 24.6 ms | 0.15 ms | 192x |
| chromosome 21, perBatch | 63 | 174 ms | 40 ms | 6.1x |
| E. coli, perBatch | 5 | | | about 1x |

(Per-update values from the timing files; the speedup uses the end-to-end
times of section 2.)


## 6. Changes to the measurement (results before them are not comparable)

Results produced before a change below measure something else and must
not be mixed with later ones. All dates 2026.

| Date, commit | Change | Effect on earlier results |
|---|---|---|
| 09.10, efa7842 | the tool's GFA gets `S` lines for new segments | before, BubbleGun skipped every link to a new segment (E. coli: 2,141 of 2,236 added edges) and analysed almost only the base graph |
| 09.10, 7a5da8e | all three BubbleGun bubble types counted; bubble lists saved | before, only the "super" type was counted (E. coli final graph: 220 of 825) |
| 09.10, edcc2df / 5c8af2c | pipelineRunner times one algorithm per call, compilation excluded | from 10.04 to 09.10 pipelineRunner did not run at all; earlier versions timed 100 runs plus compilation |
| 09.10, ecbc342 | BubbleGun timed in-process | before, every run included starting Python and importing BubbleGun (about 50 of about 70 ms on E. coli): BubbleGun appeared up to about 3x slower |
| 09.10, ac66a81 / a19887b | debug output only with `make DEBUG=1` | before, debug output (graph, orders, superbubbles after every update) ran inside the timed programs |
| 09.10, 399e0f4 | all C++ objects compiled with `-O3` | before, `superbubble.cpp`, the graph classes, the AVL tree and the AVL incremental DFS were compiled without optimisation: the algorithms were about 5x slower (e.g. static 1.54 -> 0.30 s on E. coli) |
| 10.10, abe8796 | `PeakMemory` from `VmHWM` | before, a program started by the Python runners reported the runner's peak memory when it was larger (all E. coli algorithms 13,440 KB; 27 April human files 12,032 KB) |
| 10.10, 7190235 | BubbleGun peak memory recorded | none before |
| 10.10, 22951e3 | `AlgorithmTime` = mean of 10 runs | before, one run |
| 10.10, 9482bb7 | graph defined by the walks: links no walk uses dropped, backward link traversals fixed | before, base graphs of human chromosomes contained links of no genome (chromosome 21: 979) |
| 10.10, 17895c1 | each element is added with the first genome that uses it | before, with the last one: a different update sequence (same final graph) |
| 11.10, a39591f / 12b65cc / a8fa60d | edge 0->v deleted when v gets its first real parent; incremental algorithms recompute the order with the static DFS when that happens in case 1 or 5 | before, superbubbles were missed after many updates (E. coli perEdge: static wrong on 971, incremental on 707 of 2,236; final graphs were right); incremental algorithms about 10% slower now (chromosome 21 perEdge) |

All results in `output/` produced before 10.10.2026 (including the April
human results and `scripts/ecoli_perEdge.csv`) are affected by several of
these changes and need to be rerun before they are compared.


## 7. Large graphs: HPRC chromosome 22 (experimental, branch `hprc-experiment`)

### 7.1 Data and preparation

- Graph: HPRC release 1.1 Minigraph-Cactus pangenome, GRCh38-based,
  chromosome 22 (`hprc-v1.1-mc-grch38.chroms/chr22.vg`, the standard
  clipped graph; not the `.d9` frequency-filtered or the `.full` unclipped
  version), converted with `vg convert -f`. 1,279,333 nodes, 1,780,739
  links, 54.0 Mbp; 512 walks: GRCh38 (one walk = the reference), CHM13 and
  44 HPRC samples x 2 haplotypes, split into contigs.
- Minigraph-Cactus graphs are bidirected and may contain cycles, the
  algorithms need a DAG: `Preprocessing/mcToDag.py` orients every node
  (GRCh38 nodes as on GRCh38, other nodes by the majority of the walks,
  each walk counted in the reference's direction since half of the
  contigs are reverse-complemented), drops links that would need one
  endpoint reversed (inversions) and self-loops, breaks cycles by
  removing DFS back edges only, and keeps the component of the
  reference. Everything removed is in `chromosome_22.stats.json`.
- Result (13.10.2026): 0.68% of the links removed (0.63% inversion links,
  2 cycle edges), 0.49% of the nodes (0.20% of the bases), 0.064% of all
  walk transitions broken; DAG with 1,273,074 nodes and 1,768,676 edges,
  18 sources, 36 sinks. Base graph = GRCh38 walk (740,099 nodes); the
  other 89 haplotypes add 1,028,578 edges.
- Updates as for the other datasets: perEdge = one edge at a time, in
  haplotype order, each haplotype's new edges in walk order; perBatch =
  one haplotype at a time.

### 7.2 What is measured, and how it differs from the other datasets

- **Only the best algorithm** (incMinIBST) is compared with BubbleGun.
  The comparison of our algorithms with each other is done on E. coli and
  the human chromosomes (sections 1-5); on chr22 it would only repeat it
  at a larger size, for several days of computing.
- **One run** instead of 10 (`AlgorithmTime`) / 100 (memory): a run takes
  hours, so repetitions add nothing to the precision.
- **BubbleGun perEdge is sampled.** One BubbleGun run on the full graph
  takes about 12 s (estimate), so one run after each of the ~1 million
  edges would take about 140 days. `superbubbleComparator.py --every N`
  runs BubbleGun on the graph after every N-th edge (evenly spaced over
  the whole update sequence) and estimates the total as
  `mean time of these samples x number of updates`
  (`TotalTimeEstimated` in `bubblegun_time.txt`). It also runs BubbleGun
  after the last edge of every haplotype, for the bubble lists and the
  memory, but these runs are not in the mean: they also list all bubbles
  and are slower (chr21: +13%). BubbleGun perBatch (89 updates) is not
  sampled.
- In the tables (`toolComparisonToLatex.py`) an estimated time and the
  speedup computed from it are marked with $\approx$, and the caption
  says how many updates were sampled (`~` in the CSV, which also has the
  column `Sampled updates (tool)`).
- **Validation of the sampling** (chromosome 21, where BubbleGun was run
  after all 6,863 edges): the mean of every N-th update differs from the
  mean of all updates by +0.1% (N = 50, 137 samples), +1.6% (N = 100),
  +0.1% (N = 250), +2.5% (N = 500, 13 samples), +1.8% (N = 1000, 6
  samples). Per-update time changes slowly along the sequence, so evenly
  spaced samples represent it well. For chr22 with N = 10,000 there are
  about 100 samples.
- Correctness of incMinIBST on chr22: `checkSuperbubbles.py --every 20000`
  (about 50 updates and the last one checked against the reference;
  checking all ~1 million is too slow for the Python reference).

### 7.3 Suggested wording for the paper

> To evaluate scalability we used the chromosome 22 graph of the HPRC
> release 1.1 Minigraph-Cactus pangenome (GRCh38-based; 1.28 million
> nodes, 1.78 million links, 90 haplotypes). Since our algorithms require
> a DAG, nodes were oriented by the reference and by the majority of the
> haplotype walks, links requiring an inversion (0.63% of all links) and
> two cycle-closing links were removed; 0.064% of the walk transitions
> were affected. The reference path forms the base graph (740,099 nodes)
> and the remaining 89 haplotypes add 1,028,578 edges. Running BubbleGun
> after each of these edges was infeasible (an estimated 140 days); we
> therefore measured it on the graphs after every 10,000th edge (N
> samples) and report its mean time per update multiplied by the number
> of updates. On chromosome 21, where BubbleGun could be run after every
> edge, this estimate deviates from the measured total by at most 2.5%.
> Inc-MinI-BST was run on all updates.

Fill in N (the `SampledUpdates` line of `bubblegun_time.txt`) and the
numbers from the tables; state in the table caption that the BubbleGun
time and the speedup for chr22 perEdge are estimates (the generated
caption does this).

### 7.4 Commands (on the server, branch `hprc-experiment`)

See the README section "HPRC Minigraph-Cactus Graphs".


## 8. Producing the comparison

1. `scripts/pipelineRunner.py` (preprocessing + `AlgorithmTime`),
2. `scripts/superbubbleRunner.py` (memory, superbubble lists, per-update
   counts of all algorithms),
3. `scripts/superbubbleComparator.py -t <tool>` (tool time, memory,
   bubbles),
4. `Preprocessing/toolComparisonToLatex.py -t <tool>` (tables).

Another tool (e.g. BubbleFinder) is measured with the same rules: it must
time only its own work on the given graph (no start-up), report peak
memory in a separate process, and list its bubbles with their two ends;
see `Preprocessing/superbubbleTools/toolInterface.py`.
