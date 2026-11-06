#ifndef SUPERBUBBLE_BATCH_H1
#define SUPERBUBBLE_BATCH_H1
#include<vector>
#include<ostream>
#include<utility>
#include"../graph/contextGraph.h"
#include"../DFS/contextdfs.h"
#include"superbubbleTree.h"

/*
  Superbubble computation for the batch heuristic (incBSTBatchH1): connects
  the incremental superbubble routine (SuperbubbleTree, superbubbleTree.cpp)
  to the batch DFS of incBSTBatch (one batch per genome).

  insertBatch(edges)  reads what the routine needs from the graph BEFORE the
                      batch: the groups of the batch (new vertices connected
                      by its edges, with the old vertices they touch; an edge
                      between two old vertices is a group of its own), which
                      old vertices were sinks / sources, which edges are
                      already there; then adds the batch with the batch DFS
                      (order, minI, maxO)
  updateSB()          applies the batch to the superbubbles
                      (SuperbubbleTree::updateBatch)
  computeSB()         full computation (base graph)

  prepareTime is the time of reading the groups (insertBatch, before the
  DFS); it belongs to the superbubble time of the batch.

  Superbubbles are counted and written like SuperbubbleCalculation, i.e.
  including <0 1> (the routine keeps it as its root), so the time files and
  superbubble lists are comparable with the other algorithms.
*/
class SuperbubbleCalculationBatchH1 {
public:
    SuperbubbleTree tree;
    long long duplicateEdges;     // edges already in the graph: nothing changes
    double prepareTime;           // reading the groups of the last batch (s)

public:
    // ranges: incBSTBatchH2, rescan only the ranges of the groups (updateBatchRanges)
    explicit SuperbubbleCalculationBatchH1(bool ranges = false) : duplicateEdges(0), prepareTime(0.0), pending(false), useRanges(ranges) {}
    void computeSB(DFSCONTEXT *sdfs, GraphContext *g);
    void insertBatch(DFSCONTEXT *sdfs, GraphContext *g, const std::vector<std::pair<int,int>> &edges);
    void updateSB(DFSCONTEXT *sdfs, GraphContext *g);
    void printSuperbubbles(GraphContext *g);
    void writeSuperbubbles(GraphContext *g, std::ostream &out);
    int numberOfSuperbubbles();
    // current superbubbles in internal ids, including <0 1>
    std::vector<std::pair<int,int>> superbubbles() const;

private:
    // the groups of the last batch, not yet applied by updateSB; old vertices
    // as internal ids, the new vertices of a path as input ids (their internal
    // ids exist only after the batch)
    bool pending;
    bool useRanges;
    std::vector<SuperbubbleTree::BatchGroup> groups;

    void readGroups(GraphContext *g, const std::vector<std::pair<int,int>> &edges);
    static SBView view(DFSCONTEXT *sdfs, GraphContext *g);
};
#endif
