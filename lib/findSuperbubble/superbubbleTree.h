#ifndef SUPERBUBBLE_TREE_H
#define SUPERBUBBLE_TREE_H
#include <vector>
#include <utility>

/*
  Incremental superbubble routine: the superbubbles with their nesting tree,
  updated per inserted edge by the five cases of the incremental DFS
  (see superbubbleTree.cpp for the rules of each case).

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
    };
    Stats stats;

    // full computation of all superbubbles (and the tree) for the current order
    void build(const SBView &V);
    // after inserting u->v (internal ids) and updating the order: caseNo is the
    // case of the incremental DFS; uWasSink / vWasSource as before the insertion
    void update(const SBView &V, int caseNo, int u, int v, bool uWasSink, bool vWasSource);
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
    void rescan(const SBView &V, int owner);
    void rescanOrFull(const SBView &V, int P);
};
#endif
