#ifndef INCREMENTAL_DFS_BATCH
#define INCREMENTAL_DFS_BATCH
#include "baseDFS.h"
#include "minImaxOdatastructure.h"
#include "incremental_topo.hpp"
#include "../tree/differentAVLTree.h"
#include <climits>
#include <chrono>
#include <vector>
#include <utility>

/*
  ===========================================================================
  Batch algorithm (incBSTBatch)
  ===========================================================================
  The graph grows by batches of edges (one haplotype in perBatch mode, one
  edge in perEdge mode). After each batch the order must again be a DFS
  order (reverse postorder of a DFS from the artificial source 0), because
  computeSB only works on a DFS order, and minI / maxO must be up to date.
  The superbubbles are then computed once per batch (by the program).

  The order is kept by inctopo::IncrementalTopoOrder (incremental_topo.hpp):
  its reorder fixes edges that point backwards with bounded forward searches,
  for the whole batch at once. A valid topological order is not enough,
  though: the edges that do NOT point backwards but change the DFS tree are
  handled by the same rules as in incMinIBST. Per edge u -> v of a batch
  (cases as in the per-edge algorithms, from the state at that edge):

    case 1  u new          u is placed where the DFS would visit it:
                           directly after 0 (last child of 0), or directly
                           before v if v's only parent was 0 (u takes v's
                           place; v's subtree hangs under u unchanged)
    case 2  v new          v directly after u (last child of u)
    case 3  u and v new    u, v directly before the sink 1
    case 4  u after v      the edge points backwards: batch reorder
    case 5  u before v     nothing, unless v's only parent was 0 ("adoption",
                           see below)

  Adoption (case 5, v a source): v and its subtree (a contiguous block of
  the order, from v up to the next child of 0) must hang under u. The block
  is moved directly after u (v = last child of u); then every edge from a
  vertex behind the block into the block points backwards (that vertex
  reaches part of v's old subtree first now) and is repaired by the batch
  reorder, exactly like a case 4 edge. This replaces the static DFS.
  Example: order 0 a u b v x 1, edges a->u, b->x, v->x, new edge u->v:
    move v's block after u:   0 a u v x b 1
    backward edge b->x:       reorder moves x behind b -> 0 a u v b x 1
  which is the order a fresh DFS gives (x now belongs to b).

  New vertices are not inserted one by one (each insertion shifts every
  position behind it, O(n)): step 1 only NOTES where each belongs - after or
  before an anchor vertex (u, v, 0 or 1; the anchor may itself be a new
  vertex of the batch) - and materialize() builds them all into the order in
  one pass over it, O(n + new vertices). The result is the same order as
  inserting them one at a time: several new vertices after one anchor end up
  newest first (each was inserted directly after the anchor), several before
  one anchor in the order they were noted.
  Until then a noted vertex has a virtual position directly after / before
  its anchor (vkey), exact relative to every placed vertex; comparing two
  noted vertices builds them in first.

  Steps of insertBatch():
    1. per edge: classify, add to the graph, note new vertices (cases 1-3);
       collect the new edges for the reorder; adoptions (case 5, v source)
       flush first (build in the noted vertices, reorder: the block must be
       read from a DFS order), then move the block and collect the backward
       edges
    2. flush: build in the noted vertices (one pass), then one reorder for
       all pending new edges and collected backward edges
       (IncrementalTopoOrder::insertAndRepair)
    3. copy the changed positions into topOrd / reverseTopologicalOrder
    4. minI / maxO of the vertices whose relative order changed (see
       updateMinIMaxO)

  minI is kept as in incMinIBST: every vertex c has a balanced search tree
  (minIbst[c]) of its real parents, ordered by position, and minI[c] is its
  minimum. A tree is ordered by reverseTopologicalOrder, which is only
  updated in step 3, so the vertices that change their relative order (moved
  by a reorder or an adoption) are taken out of their children's trees
  before step 3 and put back after it. maxO is kept without a tree, by
  comparisons (as in incMinIBST).
  Compile switch: with -DBATCH_MINI_SCAN (make ... MINI_SCAN=1) minI is kept
  without trees instead: a child's parents are scanned only when its
  earliest parent moved later or a parent moved earlier.
  ===========================================================================
*/
class IncrementalDFSTREEBatch : public BaseDFS, public IDatastructure {
public:
    struct IncrementalUpdateTimes {
        double insertTime = 0.0;            // step 1 without the order work: graph, classifying, noting
        double dfsTime = 0.0;               // order: building in new vertices, reorders, block moves
        double topologicalOrderTime = 0.0;  // step 3
        double minITime = 0.0;
        double maxOTime = 0.0;
        double currentinsertTime = 0.0;
        double currentdfsTime = 0.0;
        double currenttopologicalOrderTime = 0.0;
        double currentminITime = 0.0;
        double currentmaxOTime = 0.0;
    };
    struct TypeOfEdgeInsertion {
        int case1 = 0, case2 = 0, case3 = 0, case4 = 0, case5 = 0;
        int currentcase1 = 0, currentcase2 = 0, currentcase3 = 0, currentcase4 = 0, currentcase5 = 0;
    };
    struct CountOfVertices {
        int countAffectedVertices = 0;        // positions that changed
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
    TypeOfEdgeInsertion edgeTypes;
    CountOfVertices vertexCount;
    long long duplicateEdges = 0;     // edges already in the graph (ignored)
    long long adoptions = 0;          // case 5 with v a source (block moved)
    long long repairedEdges = 0;      // backward edges collected by adoptions
    long long materializations = 0;   // passes that built noted new vertices into the order

    inctopo::IncrementalTopoOrder topo;
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
    // base graph: only adds the edge to the graph (the order comes from dfs)
    void insert(GraphContext *g, int u, int v) override;
    // static DFS of the base graph; IncrementalTopoOrder starts from its order
    void dfs(GraphContext *g) override;
    // collects an edge (input ids) for the next insertBatch
    void addEdgeToBatch(int u, int v) { this->pending.push_back(std::make_pair(u, v)); }
    // applies the collected edges (steps 1-4 above)
    void insertBatch(GraphContext *g) override;
    void reset() override {}
    // full minI / maxO with the rules of SimpleDFS (base graph)
    void calculateMinI(GraphContext *g) override;
    void calculateMaxO(GraphContext *g) override;
    void printMinI() const override;
    void printMaxO() const override;

private:
    std::vector<std::pair<int,int>> pending;    // edges of the current batch (input ids)

    // ---- per batch (reset at the start of insertBatch)
    std::vector<inctopo::Edge> newEdges;        // edges for the next flush (incl. artificial ones of new vertices)
    std::vector<inctopo::Edge> backEdges;       // backward edges collected by adoptions
    std::vector<std::pair<int,int>> realEdges;  // all edges of the batch (internal ids), for minI / maxO
    std::vector<int> movedLater;                // moved by a reorder (only later relative to the others)
    std::vector<int> movedEarlier;              // moved by an adoption (a block, earlier relative to the others)
    int dirtyLo, dirtyHi;                       // positions that changed (step 3)

    // noted new vertices (step 1), built into the order by materialize()
    std::vector<std::vector<int>> afterList;    // per anchor: new vertices noted after it (in noting order)
    std::vector<std::vector<int>> beforeList;   // per anchor: new vertices noted before it (in noting order)
    std::vector<int> anchors;                   // anchors with a non-empty list (to clear them)
    std::vector<int> pendAnchor;                // per noted vertex: its anchor
    std::vector<char> pendBefore;               // per noted vertex: 1 before its anchor, 0 after
    int pendingCount = 0;

    // stamps: set membership per batch without clearing
    std::vector<unsigned> stampLater, stampEarlier, stampIn, stampOut, stampNew;
    unsigned stamp = 0;

    void fitSize(int n);
    void beginBatch();
    void notePlacement(int caseNo, int ru, int rv, bool vWasSource);
    void noteVertex(int x, int anchor, bool before);
    void materialize();
    bool isPending(int x) const { return x >= this->topo.size(); }
    long long vkey(int x) const;
    void adoptSubtree(GraphContext *g, int u, int v);
    void flushReorder();
    void markDirty(int lo, int hi);
    void copyOrder();
    void updateMinIMaxO(GraphContext *g);
#ifndef BATCH_MINI_SCAN
    void detachMovedFromTrees(GraphContext *g);   // before step 3
    void attachMovedToTrees(GraphContext *g);     // after step 3
#endif
    void minIOf(GraphContext *g, int x);
    void maxOOf(GraphContext *g, int x);
};
#endif
