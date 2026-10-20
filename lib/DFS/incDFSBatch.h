#ifndef INCREMENTAL_DFS_BATCH
#define INCREMENTAL_DFS_BATCH
#include "baseDFS.h"
#include "minImaxOdatastructure.h"
#include "incremental_topo.hpp"
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

  Steps of insertBatch():
    1. per edge: classify, add to the graph, place new vertices (cases 1-3);
       collect the new edges for the reorder; adoptions (case 5, v source)
       flush the pending reorder first (the block must be read from a DFS
       order), then move the block and collect the backward edges
    2. flush: one reorder for all pending new edges and collected backward
       edges (IncrementalTopoOrder::insertAndRepair)
    3. copy the changed positions into topOrd / reverseTopologicalOrder
    4. minI / maxO of the vertices whose relative order changed (see
       updateMinIMaxO)
  ===========================================================================
*/
class IncrementalDFSTREEBatch : public BaseDFS, public IDatastructure {
public:
    struct IncrementalUpdateTimes {
        double insertTime = 0.0;            // step 1 without reorders/moves: graph, placing
        double dfsTime = 0.0;               // reorders and block moves
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

    inctopo::IncrementalTopoOrder topo;

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

    // stamps: set membership per batch without clearing
    std::vector<unsigned> stampLater, stampEarlier, stampIn, stampOut;
    unsigned stamp = 0;

    void fitSize(int n);
    void beginBatch();
    void placeNewVertex(int caseNo, int ru, int rv, bool vWasSource);
    void adoptSubtree(GraphContext *g, int u, int v);
    void flushReorder();
    void markDirty(int lo, int hi);
    void copyOrder();
    void updateMinIMaxO(GraphContext *g);
    void minIOf(GraphContext *g, int x);
    void maxOOf(GraphContext *g, int x);
};
#endif
