# Algorithm names

Every superbubble algorithm in this repository has **one ID**. The same ID is
used everywhere the algorithm appears, so a name in a file, a command, a CSV
column or a plot always points to the same code.

| ID | Display name (plots, tables) | Topological order | minI / maxO |
|---|---|---|---|
| `static` | Static | full DFS after every update | recomputed from scratch after every update ; perEdge estimator for graphs too large to recompute after every edge (HPRC): `test/static <out> perEdge <every>` computes only after every `<every>`-th edge and after the last edge of every genome, and writes the estimated total (mean of the evenly spaced samples x updates), as `superbubbleComparator.py --every` does for the tools |
| `incRecompute` | Inc-Recompute | incremental DFS | recomputed from scratch after every update |
| `incNormalBST` | Inc-Normal-BST | incremental DFS (the same as incBST); perEdge only (in perBatch it re-keys after every edge and is far slower than every other algorithm) | both maintained in normal balanced binary search trees keyed by the position number (`lib/tree/avlTree.h`): every vertex whose position changes (cases 1 and 2: all vertices after the new one; case 4: the whole window) is taken out of its neighbours' trees and put back with its new position |
| `incBST` | Inc-BST | incremental DFS | both maintained in modified balanced binary search trees (`lib/tree/differentAVLTree.h`): they store vertex ids and compare by the current positions, so vertices that move without changing their relative order need no update |
| `incMinIBST` | Inc-MinI-BST | incremental DFS | minI maintained in an AVL tree; maxO updated directly, without a tree |
| `incRecomputeBatch` | Inc-Recompute-Batch | batch DFS as in incBSTBatch (perBatch only) | recomputed from scratch after every genome (`IncrementalDFSTREEBatch(g, true)`), as incRecompute after every edge; superbubbles once per batch |
| `incBSTBatch` | Inc-BST-Batch | batch DFS, one batch per haplotype (perBatch only; base graph: static DFS): new vertices placed as in incMinIBST in one pass (a source of the genome directly under its first parent in the genome), then one reorder for all backward edges of the batch (bounded searches, one relabel), plus one per adoption of an old source (`lib/DFS/incDFSBatch.h`); its time file has NewVertices SourcesPlaced BackwardEdges Rounds Adoptions per genome instead of the per-edge cases | updated per batch for the vertices whose relative order changed: minI in AVL trees, maxO by comparison (as in incMinIBST); superbubbles once per batch |
| `incBSTBatchH1` | Inc-BST-Batch-H1 | batch DFS as in incBSTBatch (perBatch only) | as incBSTBatch; superbubbles with the heuristic per batch (`lib/findSuperbubble/superbubbleBatchH1.cpp`, `SuperbubbleTree::updateBatch`): groups of new vertices, P = LCA on the old tree, rescan of the outermost P's, otherwise a full rebuild. Correct, but no real speedup: P = root in almost every genome |
| `incBSTBatchH2` | Inc-BST-Batch-H2 | batch DFS as in incBSTBatch (perBatch only) | as incBSTBatch; superbubbles with the heuristic per batch by ranges (`SuperbubbleTree::updateBatchRanges`): a closed group rescans only the part of P between its ends, ranges of the same P are merged, so P = root needs no full computation. Full computation only if a tail was a sink or a head a source, a dangling group is not a simple path, a closed group has a new vertex without a child or parent, or a range does not fit the tree. VERIFY=1: 0 differences on all edge lists on the laptop (not yet HPRC release 2) |

All algorithms compute the superbubbles themselves with the same code,
`SuperbubbleCalculation::computeSB` (`lib/findSuperbubble/superbubble.cpp`);
they differ only in how the topological order and minI/maxO are kept up to date.

Tables and plots list the algorithms in the order in which they were
developed (static, incRecompute, incNormalBST, incBST, incMinIBST,
incMinIBSTH1, incRecomputeBatch, incBSTBatch, incBSTBatchH1).

Which algorithms run in which mode (`BATCH_ONLY` / `PEREDGE_ONLY` in
`scripts/superbubbleRunner.py`; the parser reads the same):
- perEdge: static, incRecompute, incNormalBST, incBST, incMinIBST, incMinIBSTH1
- perBatch: static, incMinIBSTH1 (the best edge-by-edge algorithm, run per
  genome), incRecomputeBatch, incBSTBatch, incBSTBatchH1

The data structures for minI / maxO (incNormalBST, incBST, incMinIBST) are
also compared on their own (perEdge): `table_datastructures.tex` of
`Preprocessing/paperTables.py` and
`plotChromosomeComparison.py -m perEdge -a incNormalBST incBST incMinIBST -t datastructures`.

`static` is the baseline. The best edge-by-edge algorithm at the moment is
`incMinIBSTH1` (incMinIBST with the superbubble heuristic); it is the one
`scripts/pipelineRunner.py` runs by default (`BEST_ALGORITHM`). In perBatch
mode the best algorithm is `incBSTBatch` (`-a incBSTBatch`).


## Where the ID is used

Example: `incMinIBST`

| Place | Pattern | Example |
|---|---|---|
| Source file | `algorithms/<ID>.cpp` | `algorithms/incMinIBST.cpp` |
| Class in that file | `<ID with capital first letter>Superbubble` | `IncMinIBSTSuperbubble` |
| Make target | `make TARGET=<ID>` | `make TARGET=incMinIBST` |
| Makefile variables | `<ID in capitals>_SRCS/_OBJS/_EXEC` | `INCMINIBST_SRCS` |
| Binary | `test/<ID>` | `test/incMinIBST` |
| `superbubbleRunner.py` name (`-a`) | `<ID>` | `-a incMinIBST` |
| `pipelineRunner.py` (`-a`, `BEST_ALGORITHM`) | `<ID>` | `BEST_ALGORITHM = "incMinIBST"` |
| Timing file (superbubbleRunner) | `output/<dataset>/timeFile/<mode>/<ID>_time.txt` | `output/Ecoli/timeFile/perEdge/incMinIBST_time.txt` |
| Superbubble list, final graph (superbubbleRunner, normal builds) | `output/<dataset>/timeFile/<mode>/<ID>_superbubbles.txt`: `# Final graph \| total=<n>`, then `<entrance> <exit>` per line in input ids | `output/Ecoli/timeFile/perEdge/incMinIBST_superbubbles.txt` |
| Timing file (pipelineRunner) | `output/<dataset>[/<chromosome>]/incsuperbubble/<mode>/<ID>_time.txt` | `output/Ecoli/incsuperbubble/perEdge/incMinIBST_time.txt` |
| CSV columns (`parse_to_csv.py`) | `<ID><Metric>` | `incMinIBSTTotalTime`, `incMinIBSTAffectedVertices` |
| LaTeX table columns (`csvToLatex.py`) | `<ID>`, `<ID>_Speedup`, `<ID>_Imp(%)` in the data; display name in the printed header | `Inc-MinI-BST Speedup` |
| Plot legends | display name | `Inc-MinI-BST` |

`<mode>` is `perEdge` or `perbatch`. Display names are defined in
`DISPLAY_NAMES` in `Preprocessing/plotcsv.py` and `Preprocessing/csvToLatex.py`.
The chromosome comparison plots (`plotChromosomeComparison*.py`) label
`incMinIBST` as "Final" / "FINAL".


## Library classes behind each algorithm (names not changed)

The DFS / data-structure classes in `lib/DFS/` kept their original names:

| ID | DFS class | File |
|---|---|---|
| `static` | `SimpleDFS` | `lib/DFS/simpleDFS.cpp` |
| `incRecompute` | `IncrementalDFSSMM` | `lib/DFS/incrementalDFSSMM.cpp` |
| `incNormalBST` | `IncrementalDFSNormalBST` | `lib/DFS/incrementalDFSNormalBST.cpp` |
| `incBST` | `IncrementalDFSTREE` | `lib/DFS/incrementalDFSAVLTree.cpp` |
| `incMinIBST` | `IncrementalDFSTREEH4` | `lib/DFS/incrementalDFSAVLTreeH4.cpp` |
| `incBSTBatch` | `IncrementalDFSTREEBatch` | `lib/DFS/incDFSBatch.cpp` |


## Old names (before the rename)

To read older commits, notes or result files:

| ID | Old source file | Old class | Old make target | Old binary | Old runner name | Old timing file | Old CSV prefix |
|---|---|---|---|---|---|---|---|
| `static` | `superbubbledfs.cpp` | `SimpleSuperbubble` | `simpleSuperbubble` | `superbubbledfs` | `SimpleSuperbubble` | `SuperbubbleSimpleTimeOutput.txt` | `Simple` |
| `incRecompute` | `superbubbleIncdfsSMM.cpp` | `IncrementalSuperbubbleSMM` | `incsuperbubblessm` | `superbubbleIncdfsssm` | `incsuperbubblessm` | `SuperbubbleIncrementalDFSOutputSSM.txt` | `INC` |
| `incBST` | `mnlognsuperbubble.cpp` | `IncrementalSuperbubbleDifferentAVLTree` | `newincrementalSuperbubble` | `mnlognsuperbubble` | `AnotherAlgorithm` | `SuperbubbleIncrementalAVLDFSOutput.txt` | `BST` |
| `incMinIBST` | `mnlognsuperbubbleH4.cpp` | `IncrementalSuperbubbleDifferentAVLTree` | `incsuperbubbleDifferentMaxO` | `mnlognsuperbubbleH4` | `incsuperbubbleDifferentMaxO` | `SuperbubbleDifferentMaxO.txt` | `DIFFMAXO` |
| `incBSTBatch` | `superbubbleIncdfsBatch.cpp` | `IncrementalSuperbubbleBatch` | `incsuperbubblebatch` | `superbubbleIncdfsBatch` | `IncrementalSuperbubbleBatch` | `SuperbubbleIncrementalBatchDFSOutput.txt` | `INCBATCH` |

Result files written with the old names are not read by the renamed
`parse_to_csv.py`; rerun the algorithms to produce `<ID>_time.txt` files.


## Not renamed (legacy)

These are kept under their old names and are not used by the scripts:

- `algorithms/superbubbleIncdfs.cpp` (make target `incsuperbubble`): the
  original incremental DFS; commented out in `superbubbleRunner.py`
- `algorithms/superbubbledfsh1.cpp` (make target `simpleSuperbubbleh1`): H1
  heuristic version of the static algorithm
- make targets `H1Superbubble` (source `algorithms/mnlognH1.cpp` no longer
  exists) and `incdfsH1`


## External tools

External tools are not algorithms of this repository and keep their own
names: their results go to `output/<dataset>[/<chromosome>]/<tool>/<mode>/`
(e.g. `output/Ecoli/bubblegun/perEdge/bubblegun_time.txt`), where `<tool>` is
the tool's `name` in `Preprocessing/superbubbleTools/` (see
`scripts/superbubbleComparator.py`).

`Preprocessing/toolComparisonToLatex.py -t <tool>` compares the best
algorithm with any such tool (time, peak memory, and the superbubbles of
the final graph) and writes LaTeX tables to `output/toolComparison/<tool>/`.


## Build switches

- `make TARGET=<ID>`: normal build, writes the timing file, no debug output
- `make TARGET=<ID> BENCHMARK=1`: no timing file and no per-step timing (used by `pipelineRunner.py`)
- `make TARGET=<ID> DEBUG=1`: adds debug output (graph, orders, minI/maxO, superbubbles after every update); run `make clean` first
