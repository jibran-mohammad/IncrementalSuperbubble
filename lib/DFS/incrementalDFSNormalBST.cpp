#include "incrementalDFSNormalBST.h"
#include <algorithm>
#include <iostream>

/*
  Timing: TIMED(step) measures the rest of the enclosing block and adds it to
  latestTimes.<step>Time and latestTimes.current<step>Time (the same steps as
  in incBST: insert, dfs, topologicalOrder, minI, maxO). Nothing is measured
  in the BENCHMARK build.
*/
#ifndef BENCHMARK
namespace {
struct StepTimer {
    double &current, &total;
    std::chrono::high_resolution_clock::time_point start;
    StepTimer(double &c, double &t) : current(c), total(t), start(std::chrono::high_resolution_clock::now()) {}
    ~StepTimer() {
        double s = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
        current += s;
        total += s;
    }
};
}
#define TIMED(step) StepTimer stepTimer(this->latestTimes.current##step##Time, this->latestTimes.step##Time)
#define COUNT(x) x
#else
#define TIMED(step)
#define COUNT(x)
#endif


IncrementalDFSNormalBST::IncrementalDFSNormalBST(GraphContext *g) : BaseDFS(g) {
    this->cases = 0;
    g->insertEdge(0, 1);
    this->n = g->displayNumberVertices();
    growTo(this->n);
    // 0 is the first and 1 the last vertex of the order
    this->reverseTopologicalOrder[0] = 0;  this->topOrd[0] = 0;
    this->reverseTopologicalOrder[1] = 1;  this->topOrd[1] = 1;
    this->minI[0] = INT_MAX;  this->minI[1] = 0;
    this->maxO[0] = 1;        this->maxO[1] = INT_MAX;
}

void IncrementalDFSNormalBST::growTo(int size) {
    this->visit.resize(size);
    this->topOrd.resize(size);
    this->reverseTopologicalOrder.resize(size);
    this->minI.resize(size);
    this->maxO.resize(size);
    this->minIbst.resize(size);
    this->maxObst.resize(size);
    this->marked.resize(size);
    this->reachable.resize(size);
    this->postorderChanged.resize(size);
    this->intermediatepostorder.resize(size);
}


/*
  Inserts the edge u -> v (input ids) and updates the order, minI and maxO.
  The cases are those of incBST (lib/DFS/incrementalDFSAVLTree.cpp).
*/
void IncrementalDFSNormalBST::insert(GraphContext *g, int inU, int inV) {
    #ifndef BENCHMARK
    this->latestTimes.currentinsertTime = 0.0;
    this->latestTimes.currentdfsTime = 0.0;
    this->latestTimes.currenttopologicalOrderTime = 0.0;
    this->latestTimes.currentminITime = 0.0;
    this->latestTimes.currentmaxOTime = 0.0;
    this->edgeTypes.currentcase1 = this->edgeTypes.currentcase2 = this->edgeTypes.currentcase3 = 0;
    this->edgeTypes.currentcase4 = this->edgeTypes.currentcase5 = 0;
    this->vertexCount.currentcountAffectedVertices = 0;
    this->vertexCount.currentcountExploredVertices = 0;
    this->vertexCount.currentTotalAffectedEdges = 0;
    this->vertexCount.currentTotalExploredEdges = 0;
    #endif
    int newVertices;
    {
        TIMED(insert);
        int before = this->n;
        g->insertEdge(inU, inV);
        this->n = g->displayNumberVertices();
        newVertices = this->n - before;
        this->u = g->strategy->reverseMapping[inU];
        this->v = g->strategy->reverseMapping[inV];
        if (newVertices > 0) {
            growTo(this->n);
        }
        COUNT(this->vertexCount.currentTotalEdges = g->displayNumberEdges());
        COUNT(this->vertexCount.totalEdges = this->vertexCount.currentTotalEdges);
    }

    if (newVertices == 1 && g->strategy->newu) {
        this->cases = 1;
        COUNT(this->edgeTypes.case1++; this->edgeTypes.currentcase1 = 1);
        placeNewSource(g);
        g->strategy->newu = false;
    }
    else if (newVertices == 1) {
        this->cases = 2;
        COUNT(this->edgeTypes.case2++; this->edgeTypes.currentcase2 = 1);
        placeNewChild(g);
        g->strategy->newv = false;
    }
    else if (newVertices == 2) {
        this->cases = 3;
        COUNT(this->edgeTypes.case3++; this->edgeTypes.currentcase3 = 1);
        placeNewPair();
        g->strategy->newu = false;
        g->strategy->newv = false;
    }
    else if (this->reverseTopologicalOrder[u] > this->reverseTopologicalOrder[v]) {
        this->cases = 4;
        COUNT(this->edgeTypes.case4++; this->edgeTypes.currentcase4 = 1);
        reorder(g);
    }
    else {
        this->cases = 5;
        COUNT(this->edgeTypes.case5++; this->edgeTypes.currentcase5 = 1);
        forwardEdge(g);
    }
}


/*
  Case 1, u is new. As in incBST: u goes directly after 0 (position 1), or,
  if v's only parent was 0, directly before v (u takes v's place). Every
  vertex from that position on moves one position on, so it is re-keyed in
  its neighbours' trees.
*/
void IncrementalDFSNormalBST::placeNewSource(GraphContext *g) {
    int pos = g->strategy->artificialParentRemoved ? this->reverseTopologicalOrder[v] : 1;
    std::vector<int> moved = verticesAt(pos, this->n - 2);   // u is not placed yet
    removeFromTrees(g, moved);
    {
        TIMED(topologicalOrder);
        shiftOneRight(pos);
        this->reverseTopologicalOrder[u] = pos;
        this->topOrd[pos] = u;
    }
    addToTrees(g, moved);
    // u's only edge is u -> v; u has no parent but 0
    addEdgeToTrees(u, v);
    this->minI[u] = 0;
}


/*
  Case 2, v is new: v goes directly after u. Every vertex after u moves one
  position on and is re-keyed.
*/
void IncrementalDFSNormalBST::placeNewChild(GraphContext *g) {
    int pos = this->reverseTopologicalOrder[u] + 1;
    std::vector<int> moved = verticesAt(pos, this->n - 2);   // v is not placed yet
    removeFromTrees(g, moved);
    {
        TIMED(topologicalOrder);
        shiftOneRight(pos);
        this->reverseTopologicalOrder[v] = pos;
        this->topOrd[pos] = v;
    }
    addToTrees(g, moved);
    // v's only edge is u -> v; v has no child but 1
    addEdgeToTrees(u, v);
    this->maxO[v] = 1;
}


/*
  Case 3, u and v are new: u, v directly before the sink 1. Only 1 moves,
  and 1 is never stored in a tree.
*/
void IncrementalDFSNormalBST::placeNewPair() {
    {
        TIMED(topologicalOrder);
        int last = this->n - 1;
        this->reverseTopologicalOrder[1] = last;      this->topOrd[last] = 1;
        this->reverseTopologicalOrder[v] = last - 1;  this->topOrd[last - 1] = v;
        this->reverseTopologicalOrder[u] = last - 2;  this->topOrd[last - 2] = u;
    }
    addEdgeToTrees(u, v);
    this->minI[u] = 0;
    this->maxO[v] = 1;
}


/*
  Case 4, u after v (backward edge). As in incBST: partialDFS from v finds
  what must move behind u, and the window [v .. u] is relabeled. Every
  vertex of the window can get a new position, so all of them are re-keyed
  (incBST only re-inserts the vertices reached by partialDFS). The edge
  u -> v itself is in the graph already, so it goes into the trees with the
  window.
*/
void IncrementalDFSNormalBST::reorder(GraphContext *g) {
    {
        TIMED(dfs);
        partialDFS(g, v);
        COUNT(this->vertexCount.totalExploredEdges += this->vertexCount.currentTotalExploredEdges);
    }
    std::vector<int> window = verticesAt(this->reverseTopologicalOrder[v], this->reverseTopologicalOrder[u]);
    removeFromTrees(g, window);
    relabelWindow(g);
    addToTrees(g, window);
    {
        TIMED(topologicalOrder);
        reset();
    }
}


/*
  Case 5, u before v: the order stays. The edge goes into the two trees.
  If v's only parent was 0, v's subtree must hang under u: the order is
  recomputed with the static DFS (as in incBST).
*/
void IncrementalDFSNormalBST::forwardEdge(GraphContext *g) {
    if (g->strategy->artificialParentRemoved) {
        rebuildFromStaticDFS(g);
    }
    else {
        addEdgeToTrees(u, v);
    }
}


// ---------------------------------------------------------------------------
// order changes
// ---------------------------------------------------------------------------

// every vertex at a position >= from moves one position on (the new vertex
// is not in the order yet, so the last used position is n - 2)
void IncrementalDFSNormalBST::shiftOneRight(int from) {
    for (int i = this->n - 2; i >= from; i--) {
        this->topOrd[i + 1] = this->topOrd[i];
        this->reverseTopologicalOrder[this->topOrd[i + 1]] = i + 1;
    }
}

// the vertices reached from z that lie before u (as in incBST), in postorder
void IncrementalDFSNormalBST::partialDFS(GraphContext *g, int z) {
    if (this->marked[z] == 0 && this->reverseTopologicalOrder[z] < this->reverseTopologicalOrder[u]) {
        this->marked[z] = 1;
        this->reachable[z] = 1;
        COUNT(this->vertexCount.currentcountExploredVertices++);
        COUNT(this->vertexCount.currentTotalExploredEdges += g->getOutDegree(z));
        for (int w : g->displayOutgoingNeighbours(z)) {
            partialDFS(g, w);
        }
        this->m++;
        this->postorderChanged[z] = this->m;
    }
}

// the new positions of the window [v .. u] (the same as incBST's update):
// the reached vertices go behind u in postorder, the others move forward
void IncrementalDFSNormalBST::relabelWindow(GraphContext *g) {
    TIMED(topologicalOrder);
    int top = this->reverseTopologicalOrder[u];
    int shift = this->m;
    for (int i = top; i >= this->reverseTopologicalOrder[v]; i--) {
        int x = this->topOrd[i];
        COUNT(this->vertexCount.currentcountAffectedVertices++);
        COUNT(this->vertexCount.currentTotalAffectedEdges += g->getOutDegree(x));
        if (this->marked[x] == 1) {
            this->intermediatepostorder[x] = top - this->postorderChanged[x] + 1;
            this->marked[x] = 0;
        }
        else {
            this->intermediatepostorder[x] = top - shift;
            shift++;
        }
    }
    for (int x = 1; x < this->n; x++) {
        if (this->intermediatepostorder[x] != 0) {
            this->reverseTopologicalOrder[x] = this->intermediatepostorder[x];
            this->topOrd[this->reverseTopologicalOrder[x]] = x;
        }
    }
    COUNT(this->vertexCount.totalAffectedEdges += this->vertexCount.currentTotalAffectedEdges);
    COUNT(this->vertexCount.countAffectedVertices += this->vertexCount.currentcountAffectedVertices);
    COUNT(this->vertexCount.countExploredVertices += this->vertexCount.currentcountExploredVertices);
}

// the whole order with the static DFS, then every tree from the new positions
void IncrementalDFSNormalBST::rebuildFromStaticDFS(GraphContext *g) {
    {
        TIMED(dfs);
        std::fill(this->visit.begin(), this->visit.end(), 0);
        this->postorderNumber = 0;
        this->dfs(g);
    }
    {
        TIMED(minI);
        for (int c = 2; c < this->n; c++) {
            this->minIbst[c] = BalancedBinaryTree();
            for (int p : g->displayIncomingNeighbours(c)) {
                if (p != 0) this->minIbst[c].addNode(this->reverseTopologicalOrder[p], p);
            }
            updateMinI(c);
        }
    }
    {
        TIMED(maxO);
        for (int p = 2; p < this->n; p++) {
            this->maxObst[p] = BalancedBinaryTree();
            for (int c : g->displayOutgoingNeighbours(p)) {
                if (c != 1) this->maxObst[p].addNode(this->reverseTopologicalOrder[c], c);
            }
            updateMaxO(p);
        }
    }
    reset();
    this->staticRebuilds++;
}

void IncrementalDFSNormalBST::reset() {
    std::fill(this->marked.begin(), this->marked.end(), 0);
    std::fill(this->reachable.begin(), this->reachable.end(), 0);
    std::fill(this->postorderChanged.begin(), this->postorderChanged.end(), 0);
    std::fill(this->intermediatepostorder.begin(), this->intermediatepostorder.end(), 0);
    this->m = 0;
}


// ---------------------------------------------------------------------------
// the normal trees
// ---------------------------------------------------------------------------

std::vector<int> IncrementalDFSNormalBST::verticesAt(int fromPos, int toPos) const {
    std::vector<int> out;
    for (int p = fromPos; p <= toPos; p++) out.push_back(this->topOrd[p]);
    return out;
}

/*
  Before vertices move: take each of them out of its neighbours' trees, with
  its current (soon old) position as key. x is in the minI tree of each of
  its children and in the maxO tree of each of its parents.
*/
void IncrementalDFSNormalBST::removeFromTrees(GraphContext *g, const std::vector<int> &moved) {
    {
        TIMED(minI);
        for (int x : moved) {
            if (x == 0 || x == 1) continue;
            for (int c : g->displayOutgoingNeighbours(x)) {
                if (c != 1) this->minIbst[c].deleteNode(this->reverseTopologicalOrder[x]);
            }
        }
    }
    {
        TIMED(maxO);
        for (int x : moved) {
            if (x == 0 || x == 1) continue;
            for (int p : g->displayIncomingNeighbours(x)) {
                if (p != 0) this->maxObst[p].deleteNode(this->reverseTopologicalOrder[x]);
            }
        }
    }
}

// after the move: put them back with their new positions, and update minI /
// maxO of every vertex whose tree changed
void IncrementalDFSNormalBST::addToTrees(GraphContext *g, const std::vector<int> &moved) {
    {
        TIMED(minI);
        for (int x : moved) {
            if (x == 0 || x == 1) continue;
            for (int c : g->displayOutgoingNeighbours(x)) {
                if (c == 1) continue;
                this->minIbst[c].addNode(this->reverseTopologicalOrder[x], x);
                updateMinI(c);
            }
        }
    }
    {
        TIMED(maxO);
        for (int x : moved) {
            if (x == 0 || x == 1) continue;
            for (int p : g->displayIncomingNeighbours(x)) {
                if (p == 0) continue;
                this->maxObst[p].addNode(this->reverseTopologicalOrder[x], x);
                updateMaxO(p);
            }
        }
    }
}

// the edge a -> b: a into b's minI tree, b into a's maxO tree
void IncrementalDFSNormalBST::addEdgeToTrees(int a, int b) {
    {
        TIMED(minI);
        this->minIbst[b].addNode(this->reverseTopologicalOrder[a], a);
        updateMinI(b);
    }
    {
        TIMED(maxO);
        this->maxObst[a].addNode(this->reverseTopologicalOrder[b], b);
        updateMaxO(a);
    }
}

// minI[c]: the parent with the smallest position, 0 if c has no real parent
void IncrementalDFSNormalBST::updateMinI(int c) {
    this->minI[c] = this->minIbst[c].avlTree.empty() ? 0 : this->minIbst[c].findMin().index;
}

// maxO[p]: the child with the largest position, 1 if p has no real child
void IncrementalDFSNormalBST::updateMaxO(int p) {
    this->maxO[p] = this->maxObst[p].avlTree.empty() ? 1 : this->maxObst[p].findMax().index;
}


void IncrementalDFSNormalBST::printMinI() const {
    std::cout << "MinI of the graph is" << std::endl;
    for (int i = 0; i < this->n; i++) std::cout << i << " -> " << this->minI[i] << std::endl;
}

void IncrementalDFSNormalBST::printMaxO() const {
    std::cout << "MaxO of the graph is " << std::endl;
    for (int i = 0; i < this->n; i++) std::cout << i << " -> " << this->maxO[i] << std::endl;
}
