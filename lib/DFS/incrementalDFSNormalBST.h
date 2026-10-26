#ifndef INCREMENTAL_DFS_NORMAL_BST
#define INCREMENTAL_DFS_NORMAL_BST
#include "baseDFS.h"
#include "minImaxOdatastructure.h"
#include "../tree/avlTree.h"
#include <climits>
#include <chrono>
#include <vector>

/*
  incNormalBST: the incremental DFS of incBST (same five cases, same
  partialDFS and relabeling), but minI / maxO are kept in NORMAL binary
  search trees.

    minIbst[c]  the real parents of c, keyed by their position in the
                topological order  -> minI[c] = parent with the smallest key
    maxObst[p]  the real children of p, keyed by their position
                -> maxO[p] = child with the largest key

  The difference to incBST: incBST's modified tree (DifferentAVLTree) stores
  vertex ids and reads their CURRENT positions when it compares, so a vertex
  can move without touching any tree as long as the relative order stays
  the same. A normal tree stores the position NUMBER as key; when a vertex
  moves, that key is out of date. So every vertex whose position changes is
  taken out of its neighbours' trees before the change (with the old
  position) and put back after it (with the new position):

    case 1 / 2   every vertex after the new one moves one position on
    case 4       every vertex of the window [v .. u] gets a new position
    case 3, 5    no existing vertex moves (only the sink 1, never stored)

  Vertices 0 and 1 are never stored in a tree: minI = 0 means "only the
  source 0 is a parent", maxO = 1 means "only the sink 1 is a child".
*/
class IncrementalDFSNormalBST : public BaseDFS, public IDatastructure {
public:
    struct IncrementalUpdateTimes {
        double insertTime = 0.0;
        double dfsTime = 0.0;
        double topologicalOrderTime = 0.0;
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
        int countAffectedVertices = 0;        // positions in the reorder window (case 4)
        int countExploredVertices = 0;        // vertices reached by partialDFS
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

    std::vector<BalancedBinaryTree> minIbst;   // parents of a vertex, key = position
    std::vector<BalancedBinaryTree> maxObst;   // children of a vertex, key = position
    int staticRebuilds = 0;                    // case 5 with v's only parent 0

public:
    IncrementalDFSNormalBST(GraphContext *g);
    void insert(GraphContext *g, int u, int v) override;
    void reset() override;
    void calculateMinI(GraphContext *g) override {}   // kept up to date by insert
    void calculateMaxO(GraphContext *g) override {}
    void printMinI() const override;
    void printMaxO() const override;

private:
    int u, v;                          // the inserted edge (internal ids)
    int m = 0;                         // postorder counter of partialDFS
    std::vector<int> marked, reachable, postorderChanged, intermediatepostorder;

    void growTo(int n);
    // the cases
    void placeNewSource(GraphContext *g);        // case 1
    void placeNewChild(GraphContext *g);         // case 2
    void placeNewPair();                         // case 3
    void reorder(GraphContext *g);               // case 4
    void forwardEdge(GraphContext *g);           // case 5
    // order changes (as in incBST)
    void shiftOneRight(int from);                // positions >= from move one on
    void partialDFS(GraphContext *g, int z);
    void relabelWindow(GraphContext *g);
    void rebuildFromStaticDFS(GraphContext *g);
    // the normal trees
    std::vector<int> verticesAt(int fromPos, int toPos) const;
    void removeFromTrees(GraphContext *g, const std::vector<int> &moved);
    void addToTrees(GraphContext *g, const std::vector<int> &moved);
    void addEdgeToTrees(int a, int b);           // a -> b, both with current positions
    void updateMinI(int c);
    void updateMaxO(int p);
};
#endif
