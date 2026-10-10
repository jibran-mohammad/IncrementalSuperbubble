# How the pangenome graphs and their updates are prepared

This document describes, from the code, how a dataset goes from genome
FASTA files to the sequence of graphs that the superbubble algorithms (and
external tools such as BubbleGun) process. It states the rules that define
the graphs, the checks that enforce them, the numbers for the datasets,
and what must change for a later version of the algorithm (e.g. one that
keeps orientations).

The steps are run by `scripts/graphRunner.py` (all of them) and
`scripts/pipelineRunner.py` (steps 2-6, starting from an existing graph).
Files are in `input/<dataset>/Graphs/`.

| Step | Script | Writes |
|---|---|---|
| 1 | `Preprocessing/graphgeneration.py` | `<ref>_mg.gfa` (graph without walks), `<ref>.gfa` (graph with one walk per genome) |
| 2 | `Preprocessing/normalizeOrientations.py` | `<ref>_pos.gfa` (all links and walk steps forward) |
| - | `Preprocessing/DAGCheck.py` (graphRunner only) | nothing (checks `<ref>.gfa`) |
| 3 | `Preprocessing/haplotypeRemoval.py` | `reduced_graph.gfa` (base graph), `unique_elements.txt` (what each genome adds), `haplotype_removal_report.txt` (counts) |
| 4 | `Preprocessing/haplotypeSorter.py` | `ordered_haplotypes.txt` (order of the updates) |
| 5 | `Preprocessing/gfaToEdgeListConverter.py` | `full_graph_edgelist.txt` |
| 6 | `Preprocessing/edgeListProcessor.py` | `numeric_edgelist.txt` (input of the algorithms) |


## 1. Graph and walks (`graphgeneration.py`)

1. **Order of the genomes.** The Mash distance of every genome to the
   reference is computed (`mash dist`); genomes are used in increasing
   distance (closest first).
2. **Graph.** `minigraph -cxggs` (or `minichain -cxggs`, option
   `--graph-backend`) builds the graph from the reference and adds the
   genomes in that order: each genome is aligned to the graph as it is at
   that moment and its variants are inserted as new segments and links.
   The graph has no walks; every segment and link carries the rank of the
   genome that introduced it (`SR` tag, 0 = reference). It is saved as
   `<ref>_mg.gfa`.
3. **Walks.** Every genome is aligned again, to the finished graph, with
   `minichain -cx asm --vc`. From each genome's alignments **one record is
   kept: the best primary record** (primary first, then highest mapping
   quality); its path becomes the genome's walk (one `W` line per genome).
   The reference walk comes from `vg convert` if vg runs, otherwise from
   aligning the reference the same way. The alignment files are deleted
   afterwards.
4. **Segments on no walk** are removed, and so are the links attached to
   them. Links are not removed otherwise (no link whose two segments are
   kept was removed: checked on chromosome 21).

The graph is therefore **defined by the walks**: what is analysed is the
union of one walk per genome, not the whole construction graph.


## 2. Orientation (`normalizeOrientations.py`)

The superbubble algorithms work on a directed graph in which every walk is
read forward. A segment visited backwards (`<sN` in a walk, or `sN -` in a
link) is replaced by a copy `sN'` that is visited forwards; links and walks
are rewritten accordingly, so every link is `+ -> +` and every walk step is
`>`.

A link can be traversed by a walk in either direction. Written once in
its stored direction (`a -> b'` for `a + b -`), the opposite traversal
(`b -> a'`) would have no link. **For every walk step without a link, a
forward link is added** (it is the same link of the original graph, read
backwards; all links have overlap `0M`). Chromosome 21: 11 such links;
E. coli: none.


## 3. Base graph and updates (`haplotypeRemoval.py`)

Rules (each element = a segment or a link):

1. **An element that no walk uses is dropped.** Without this rule such
   elements were never given to a genome and ended up in the base graph
   (chromosome 21: 979 links, 125 of them pointing at segments that are
   not in the base graph).
2. **Base graph = exactly the segments and links of the reference walk.**
3. **Every other element is given to one genome: the first one, in
   walk (Mash) order, whose walk uses it.** Genomes are removed from the
   full graph one at a time, starting with the last one: each element has
   a use count over all walks, removing a genome decrements the counts of
   the elements its walk uses, and an element whose count reaches 0 (no
   remaining walk uses it) belongs to the genome just removed. The
   reference is never removed, so its elements stay in the base graph.

Checks (the script stops with a message if one fails):

- one walk per genome (walks are stored by genome name);
- all links `+ -> +` and all walk steps forward (step 2 has run);
- the graph of all walks is acyclic (the algorithms need a DAG); this
  checks the graph that is actually used, after step 2, whereas
  `DAGCheck.py` checks `<ref>.gfa`;
- the written base graph equals the reference walk, and every link in it
  connects two segments it defines;
- adding the genomes in order, every genome's walk is complete in the
  graph right after that genome is added.

The counts are written to `haplotype_removal_report.txt`.

**Consequence of rule 3: the intermediate graphs.** Updates are added in
walk (Mash) order. Because every element comes with the first genome that
uses it, the graph after adding genomes 1..k is **exactly the union of the
walks of the reference and genomes 1..k**, for every k (checked, see
above). The final graph contains every walk.

An earlier version removed genomes starting with the *first* one, which
gives every element to the *last* genome that uses it. The intermediate
graphs were then not unions of whole walks: a genome's walk was complete
right after adding it for only 1 of 5 E. coli genomes (min 73.9 %) and
1 of 63 chromosome 21 genomes (min 63.2 %). The final graph and the set of
added edges are the same under both rules; only which update an edge
belongs to differs (E. coli: 607 of 2,236 edges, chromosome 21: 4,992 of
6,863). Results computed before this change are per-update results of the
other rule.


## 4. Order of the updates (`haplotypeSorter.py`)

Genomes are put in three groups, keeping their order within a group:
genomes that add edges, genomes that add only segments, genomes that add
nothing. In all datasets so far every genome adds edges, so the update
order is the walk (Mash) order.


## 5. Edge list (`gfaToEdgeListConverter.py`, `edgeListProcessor.py`)

`numeric_edgelist.txt`: a header `<base edges> <base nodes>`, the base
edges, then for every genome (in update order) a `# HAPLOTYPE: <name>` line,
its number of edges and its edges. Segment `sX` is numbered `X+1`. The
algorithms add an artificial source (0) and sink (1) and connect sources
and sinks to them.

- **perEdge mode:** one update = one edge.
- **perBatch mode:** one update = all edges of one genome.


## 6. What the algorithms and the external tools see

Both process the same graphs:

- the C++ algorithms read `numeric_edgelist.txt`;
- `scripts/superbubbleComparator.py` gives the tool, after every update,
  a GFA with the base graph (`reduced_graph.gfa`), the edges added so far
  (`unique_elements.txt`, in `ordered_haplotypes.txt` order) and an `S` line
  for every segment those edges introduce. The tool stops the run if it
  reports a warning (e.g. a link to an undefined segment).

When comparing superbubbles (`Preprocessing/toolComparisonToLatex.py`),
pairs are unordered GFA segment names; superbubbles with the artificial
source/sink and single-edge superbubbles are left out on both sides.


## 7. Datasets

| | E. coli | chromosome 21 |
|---|---|---|
| genomes (incl. reference) | 6 | 64 |
| graph without walks (`_mg.gfa`): segments / links | 2,651 / 3,645 | 23,869 / 34,021 |
| graph with walks (`.gfa`): segments / links | 2,644 / 3,631 | 14,490 / 18,556 |
| after orientation (`_pos.gfa`): segments / links | 2,644 / 3,631 | 14,629 / 18,567 |
| dropped, used by no walk: segments / links | 0 / 1 | 139 / 979 |
| graph of all walks (used): segments / links | 2,644 / 3,630 | 14,490 / 17,588 |
| base graph = reference walk: segments / links | 1,395 / 1,394 | 10,726 / 10,725 |
| added by genomes: segments / links (= perEdge updates) | 1,249 / 2,236 | 3,764 / 6,863 |

(Chromosome 21: of the 18,556 links of the graph with walks, 968 are used
by no walk in either direction; 18,556 - 968 = 17,588.)

Effect of the rules on chromosome 21 (perBatch, final graph): before
them, the base graph contained 979 links of no genome and BubbleGun had
to skip 125 links to undefined segments; the algorithms reported 9,020
(static) and 9,030 (incremental) superbubbles. With them, BubbleGun runs
without warnings, all four algorithms report 9,999 superbubbles, and
Inc-MinI-BST and BubbleGun find the same 2,573 superbubbles (single-edge
ones and the artificial source/sink excluded; see section 6).


## 8. Known limitations

- **One alignment record per genome.** Each genome's walk is its best
  single alignment record (section 1.3). If a genome aligns in several
  records (large rearrangements, assembly gaps), the rest of it is not in
  its walk, and elements only those parts use are not in the analysed
  graph. The alignment files are deleted, so this cannot be measured for
  existing graphs. On chromosome 21 the non-reference walks are 65-89 %
  (median 75 %) of the reference walk's length; part of this can also be
  real (e.g. the acrocentric short arm, often missing in assemblies).
- **Construction graph vs. analysed graph.** Segments and links of the
  construction graph that no walk uses are not analysed (chromosome 21:
  9,379 segments removed in step 1, 979 links in step 3).


## 9. For the next version (e.g. keeping orientations)

The checks in `haplotypeRemoval.py` encode the assumptions of the current
algorithms. A version that works on another kind of graph has to change
these points on purpose (the checks make the pipeline stop otherwise):

- **Orientation.** If the algorithm keeps orientations, skip
  `normalizeOrientations.py`; then links must be matched in both
  directions (a walk can use `a + b -` as `<b <a`) when counting uses and
  when giving links to genomes, and the "all forward" check must go.
- **Cycles.** If the algorithm accepts cycles, remove the acyclicity check
  (orientation-preserving graphs of chromosomes often have cycles).
- **Several walks per genome.** To use all alignment records of a genome,
  `graphgeneration.py` must write one walk per record (and keep the
  alignment files), and `haplotypeRemoval.py` must store walks per
  (genome, piece) instead of per genome; the one-walk-per-genome check
  then becomes a per-piece check.
- **Order of elements.** Keep giving each element to the first genome
  that uses it (section 3), so that intermediate graphs are unions of
  whole walks; the incremental check enforces it.
- Keep the other rules and checks: base graph = reference walk, every
  link in the base graph has both segments, every walk step is a link,
  elements used by no walk are dropped, every walk is complete when its
  genome is added.

To check a dataset: run the preprocessing (e.g. `scripts/pipelineRunner.py`)
and read `input/<dataset>/Graphs/haplotype_removal_report.txt`; the BubbleGun
comparison (`scripts/superbubbleComparator.py`) additionally fails on any
link to an undefined segment.
