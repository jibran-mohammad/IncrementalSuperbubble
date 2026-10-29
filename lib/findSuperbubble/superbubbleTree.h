#ifndef SUPERBUBBLE_TREE_H
#define SUPERBUBBLE_TREE_H
#include <vector>
#include <utility>

/*
  Incremental superbubble routine: the superbubbles with their nesting tree,
  updated per inserted edge by the five cases of the incremental DFS
  (see superbubbleTree.cpp for the rules of each case), or per batch of
  edges (one genome, incBSTBatchH1) by groups of the batch (updateBatch).

  SBView is a read-only view of the current DFS order (no copy):
    n        number of vertices (0 = artificial source, 1 = artificial sink)
    pos[x]   position of x in the topological order (reverseTopologicalOrder)
    at[p]    vertex at position p (topOrd)
    minI[x]  earliest parent of x, maxO[x] latest child of x (vertex ids)

  Nodes are kept in a pool (ns/nt/npar/nalive); node ROOT = <0 1>. Killed
  nodes stay in the pool with their parent pointers, findAlive() skips them
  (with path compression). Per vertex:
    byEntrance[x]  node with entrance x (-1: none)
    byExit[x]      node with exit x (-1: none)
    encl[x]        a node containing x strictly inside (possibly dead; the
                   innermost alive one is findAlive(encl[x]))
*/
struct SBView {
    int n;
    const int *pos;
    const int *at;
    const int *minI;
    const int *maxO;
};

class SuperbubbleTree {
public:
    enum { ROOT = 0 };

    struct Stats {
        long long added = 0, killed = 0;
        long long rescans = 0, rescanVertices = 0;
        long long fullRecomputes = 0, pathOnly = 0;
        // why full: u was a sink or v a source (4b, 5c, 5d) / P is the root (5a)
        long long fullSinkSource = 0, fullRootP = 0;
        // batches (updateBatch): batches without a scan; why full: a dangling
        // group that is not a simple path or that destroys a region to rescan /
        // a region that is no longer a superbubble after the batch
        long long batches = 0, batchesNoScan = 0;
        long long fullDangling = 0, fullNotClosed = 0;
    };

    /*
      A group of a batch (built by SuperbubbleCalculationBatchH1 from the
      graph BEFORE the batch): the new vertices that the batch's edges connect,
      with the old vertices they touch, or one edge between two old vertices.
        tails / heads    old vertices with an edge into / out of the group
        tailWasSink      a tail had no real child before the batch
        headWasSource    a head had no real parent before the batch
      Dangling group (no tails or no heads):
        isPath           its new vertices form one simple path, in 'path'
        tailPair         (tail, path front) becomes a superbubble: the tail was
                         a sink and gets only this one new child
        headPair         (path back, head) likewise
    */
    struct BatchGroup {
        std::vector<int> tails, heads;
        bool tailWasSink = false, headWasSource = false;
        bool isPath = false, tailPair = false, headPair = false;
        std::vector<int> path;
    };
    Stats stats;

    // full computation of all superbubbles (and the tree) for the current order
    void build(const SBView &V);
    // after inserting u->v (internal ids) and updating the order: caseNo is the
    // case of the incremental DFS; uWasSink / vWasSource as before the insertion
    void update(const SBView &V, int caseNo, int u, int v, bool uWasSink, bool vWasSource);
    // after a whole batch: V is the order after the batch, the tree still the
    // one before it; the groups as above (internal ids)
    void updateBatch(const SBView &V, const std::vector<BatchGroup> &groups);
    // current superbubbles (without the root <0 1>), as (entrance, exit)
    std::vector<std::pair<int,int>> superbubbles() const;
    int count() const { return aliveCount; }
    bool checkAgainstFullScan(const SBView &V) const;

private:
    std::vector<int> ns, nt, npar, nmark;
    std::vector<char> nalive;
    std::vector<int> byEntrance, byExit, encl;
    int aliveCount = 0;
    int stamp = 0;

    std::vector<std::pair<int,int>> scanRange(const SBView &V, int L, int R) const;
    void ensureSize(int n);
    int newNode(int s, int t, int parent);
    void kill(int id);
    int findAlive(int id);
    int innerOut(int x);
    int innerIn(int x);
    int lca(int a, int b);
    void killPath(int from, int to);
    bool rescan(const SBView &V, int owner);
    void rescanOrFull(const SBView &V, int P);
    void fullBuild(const SBView &V, long long &reason);
    void addLeaf(int s, int t);
};
#endif
