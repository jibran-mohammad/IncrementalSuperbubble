#ifndef SUPERBUBBLE_H1
#define SUPERBUBBLE_H1
#include<vector>
#include<ostream>
#include<utility>
#include"../graph/contextGraph.h"
#include"../DFS/contextdfs.h"
#include"superbubbleTree.h"

/*
  Superbubble computation for the heuristic (incMinIBSTH1): connects the
  incremental superbubble routine (SuperbubbleTree, superbubbleTree.cpp) to
  the incremental DFS of incMinIBST.

  insertE(u,v)  reads what the routine needs from the graph BEFORE the edge
                is added (is u a sink, is v a source, is the edge already
                there), then adds the edge with the incremental DFS (order,
                minI, maxO) and remembers its case
  updateSB()    applies that edge to the superbubbles (SuperbubbleTree::update)
  computeSB()   full computation (base graph, perBatch)

  Superbubbles are counted and written like SuperbubbleCalculation, i.e.
  including <0 1> (the routine keeps it as its root), so the time files and
  superbubble lists are comparable with the other algorithms.
*/
class SuperbubbleCalculationH1 {
public:
    SuperbubbleTree tree;
    long long duplicateEdges;     // edges already in the graph: nothing changes

public:
    SuperbubbleCalculationH1() : duplicateEdges(0), pending(false) {}
    void computeSB(DFSCONTEXT *sdfs, GraphContext *g);
    void insertE(DFSCONTEXT *sdfs, GraphContext *g, int u, int v);
    void updateSB(DFSCONTEXT *sdfs, GraphContext *g);
    void printSuperbubbles(GraphContext *g);
    void writeSuperbubbles(GraphContext *g, std::ostream &out);
    int numberOfSuperbubbles();
    // current superbubbles in internal ids, including <0 1>
    std::vector<std::pair<int,int>> superbubbles() const;

private:
    // the last inserted edge, not yet applied by updateSB
    bool pending;
    int pU, pV, pCase;
    bool pUWasSink, pVWasSource, pDuplicate;

    static SBView view(DFSCONTEXT *sdfs, GraphContext *g);
};
#endif
