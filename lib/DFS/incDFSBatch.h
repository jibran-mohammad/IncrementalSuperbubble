#ifndef INCREMENTAL_DFS_BATCH
#define INCREMENTAL_DFS_BATCH
#include "baseDFS.h"
#include "minImaxOdatastructure.h"
#include "../tree/differentAVLTree.h"
#include <climits>
#include <chrono>
#include <list>
#include <vector>
#include <utility>

/*
  ===========================================================================
  Batch algorithm (incBSTBatch)
  ===========================================================================
  The graph grows by batches of edges, one batch per genome. After each
  batch the order must again be a DFS order (reverse postorder of a DFS from
  the artificial source 0), because computeSB only works on a DFS order, and
  minI / maxO must be up to date. The superbubbles are then computed once per
  batch (by the program).

  Base graph: insert() adds every edge to the graph, dfs() is the static DFS
  (as in the static algorithm), calculateMinI / calculateMaxO build minI and
  maxO.

  A batch (SuperbubbleCalculation::insertBatch -> DFSCONTEXT::updateEdgeBatch
  -> insertBatch(g, edges)), the edges of one genome in their order:

    1. classify: one pass over the edges, without changing the graph. A
       vertex is new at an edge if it is neither in the graph nor seen at an
       earlier edge of the batch; v's only parent is 0 if it has no real
       parent in the graph and got none at an earlier edge of the batch.
         E_new  edges with a new vertex (cases 1-3), each new vertex gets
                its place next to an anchor (rules of incMinIBST):
                  case 1  u new    u directly after 0 (last child of 0), or
                                   directly before v if v's only parent was 0
                                   (u takes v's place, v's subtree hangs
                                   under u unchanged)
                  case 2  v new    v directly after u (last child of u)
                  case 3  both     u directly before the sink 1, v after u
                A new vertex that is a source when it appears (u of case 1
                or 3) gets its place only at the end of the pass: if a later
                edge x -> u of the genome gives it its first parent, u goes
                directly after x (where an adoption or a reorder would put
                it; the edge counts to E_new), otherwise to its default place
         E_old  edges between vertices that exist at that point
       Duplicate edges are dropped. New vertices get the internal ids the
       graph will give them (consecutive, in order of first appearance).
    2. new vertices: E_new goes into the graph; materialize() builds all new
       vertices into topOrd / reverseTopologicalOrder in one pass over the
       order, O(n + new vertices). Every relative order of the old vertices
       is kept.
    3. E_old, in order (each edge goes into the graph first):
         u after v (backward)   collected for the reorder (case 4)
         u before v             nothing (case 5), unless v's only parent was
                                0: adoption (below)
       repairRound(): the batch reorder for all collected backward edges at
       once (the IncrementalDFS / Update of the batch DFS algorithm): the
       window of positions from the earliest head to the latest tail is
       scanned from right to left; for every tail (in its edges' order) a
       search from the head, bounded by the tail's position, collects the
       vertices that must move behind the tail, in postorder; each vertex is
       emitted once and the window is relabeled once.
       E_new edges that point backwards after step 2 (a source placed under
       its parent with an edge to an old vertex before that parent) go into
       the first round.
       Adoption (v's only parent was 0 - with the rule above only old
       sources, rare - v's subtree must hang under u): the
       collected backward edges are repaired first (the block must be read
       from a DFS order), then v's block (v up to the next child of 0) moves
       directly after u; edges from vertices now behind the block into it
       point backwards and are collected for the next round.
       At the end, one round for the rest. Without adoptions a batch has
       exactly one round.
    4. minI / maxO of the vertices whose relative order changed (see
       updateMinIMaxO); minI in trees as in incMinIBST, maxO by comparisons.

  The minI trees compare by reverseTopologicalOrder: a vertex is taken out
  of its children's trees before its relative order changes (searched in a
  round, or in an adopted block) and put back after step 3. Building in new
  vertices keeps every relative order.
  Compile switch: with -DBATCH_MINI_SCAN (make ... MINI_SCAN=1) minI is kept
  without trees: a child's parents are scanned only when its earliest
  parent moved later or a parent moved earlier.
  ===========================================================================
*/
class IncrementalDFSTREEBatch : public BaseDFS, public IDatastructure {
public:
    struct IncrementalUpdateTimes {
        double insertTime = 0.0;            // step 1 and the graph
        double dfsTime = 0.0;               // searches of the reorder, block moves
        double topologicalOrderTime = 0.0;  // building in new vertices, relabeling
        double minITime = 0.0;
        double maxOTime = 0.0;
        double currentinsertTime = 0.0;
        double currentdfsTime = 0.0;
        double currenttopologicalOrderTime = 0.0;
        double currentminITime = 0.0;
        double currentmaxOTime = 0.0;
    };
    // what a batch did (the per-edge cases 1-5 of the other algorithms do not
    // apply: the edges of a batch are handled together); current* = this batch
    struct BatchCounts {
        long long newVertices = 0;        // new vertices, placed in one pass
        long long sourcesPlaced = 0;      // sources of the genome placed directly under their first parent
        long long backwardEdges = 0;      // edges repaired by the reorder rounds
        long long rounds = 0;             // reorder rounds (1 per batch with backward edges + 1 per adoption)
        long long adoptions = 0;          // old sources that got their first parent (block moved under it)
        long long currentNewVertices = 0, currentSourcesPlaced = 0, currentBackwardEdges = 0;
        long long currentRounds = 0, currentAdoptions = 0;
    };
    struct CountOfVertices {
        int countAffectedVertices = 0;        // positions in the reorder windows and moved blocks (affected edges: their out-degrees)
        int countExploredVertices = 0;        // vertices reached by the reorder's searches
        int currentcountAffectedVertices = 0;
        int currentcountExploredVertices = 0;
        int totalEdges = 0;
        int currentTotalEdges = 0;
        int totalExploredEdges = 0;
        int totalAffectedEdges = 0;
        int currentTotalExploredEdges = 0;
        int currentTotalAffectedEdges = 0;
    };
    IncrementalUpdateTimes latestTimes;
    BatchCounts counts;
    CountOfVertices vertexCount;
    long long duplicateEdges = 0;     // edges already in the graph (ignored)
    long long repairedEdges = 0;      // backward edges collected by adoptions

#ifndef BATCH_MINI_SCAN
    // minIbst[c]: the real parents of c, ordered by position (minI[c] = minimum)
    std::vector<DifferentBalancedBinaryTree> minIbst;
#endif

public:
    IncrementalDFSTREEBatch(GraphContext *g) : BaseDFS(g) {
        cases = 0;
        g->insertEdge(0, 1);
        n = g->displayNumberVertices();
        this->fitSize(n);
    }
    // base graph: the edge goes into the graph (the order comes from dfs)
    void insert(GraphContext *g, int u, int v) override;
    // static DFS of the base graph (as in the static algorithm)
    void dfs(GraphContext *g) override;
    // a batch: the edges of one genome (input ids), steps 1-4 above
    void insertBatch(GraphContext *g, const std::vector<std::pair<int,int>> &edges) override;
    void reset() override {}
    // full minI / maxO with the rules of SimpleDFS (base graph)
    void calculateMinI(GraphContext *g) override;
    void calculateMaxO(GraphContext *g) override;
    void printMinI() const override;
    void printMaxO() const override;

private:
    struct OldEdge { int x, y, u, v; bool vOnly0; };   // input ids, internal ids

    // ---- per batch
    std::vector<std::pair<int,int>> newInputEdges;  // E_new (input ids), in order
    std::vector<OldEdge> oldEdges;              // E_old, in order
    std::vector<std::pair<int,int>> pending;    // backward edges for the next round
    std::vector<std::pair<int,int>> realEdges;  // all edges of the batch (internal ids), for minI / maxO
    std::vector<int> movedLater;                // searched in a round (only later relative to the others)
    std::vector<int> movedEarlier;              // in an adopted block (earlier relative to the others)
    std::vector<int> detached;                  // taken out of their children's minI trees

    // new vertices, built into the order by materialize()
    std::vector<std::vector<int>> afterList;    // per anchor: new vertices noted after it
    std::vector<std::vector<int>> beforeList;   // per anchor: new vertices noted before it
    std::vector<int> anchors;                   // anchors with a non-empty list (to clear them)

    // stamps: set membership per batch without clearing
    std::vector<unsigned> stampLater, stampEarlier, stampIn, stampOut, stampNew, stampParent, stampDetached;
    unsigned stamp = 0;
    // per round: searched vertices, tails and their edges (in order)
    std::vector<unsigned> markEpoch, tailEpoch;
    std::vector<int> tailFirst, tailLast, edgeNext;
    std::vector<int> out;                       // the window, emitted from right to left
    unsigned epoch = 0;

    void fitSize(int n);
    void beginBatch();
    void noteVertex(int x, int anchor, bool before);
    void materialize(int nBefore, int nAfter);
    void adoptSubtree(GraphContext *g, int u, int v);
    void repairRound(GraphContext *g);
    long long search(GraphContext *g, int v, int tail, int bound);
    void detachFromTrees(GraphContext *g, int x);
    void updateMinIMaxO(GraphContext *g);
#ifndef BATCH_MINI_SCAN
    void attachToTrees(GraphContext *g);
#endif
    void minIOf(GraphContext *g, int x);
    void maxOOf(GraphContext *g, int x);
};
#endif
