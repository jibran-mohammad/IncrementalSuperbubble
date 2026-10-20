#include"incDFSBatch.h"
#include<algorithm>
#include<cstdlib>

/*
  Batch algorithm: see the description in incDFSBatch.h. Order of this file:
    base graph          insert, dfs, calculateMinI / calculateMaxO
    a batch             insertBatch (steps 1-4) and its helpers in step order:
                        beginBatch, placeNewVertex, adoptSubtree,
                        flushReorder, copyOrder, updateMinIMaxO
    minI / maxO of one vertex (rules of SimpleDFS), printing
*/

#ifndef BENCHMARK
#define BATCH_TIMER_START auto batchTimerStart = std::chrono::high_resolution_clock::now();
#define BATCH_TIMER_ADD(var) var += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - batchTimerStart).count();
#else
#define BATCH_TIMER_START
#define BATCH_TIMER_ADD(var)
#endif

    // the per-vertex arrays grow with the graph
    void IncrementalDFSTREEBatch::fitSize(int n) {
        if ((int)this->topOrd.size() < n) {
            this->visit.resize(n, 0);
            this->topOrd.resize(n, 0);
            this->reverseTopologicalOrder.resize(n, 0);
            this->oldReverseTopologicalOrder.resize(n, 0);
            this->minI.resize(n, 0);
            this->maxO.resize(n, 0);
            this->stampLater.resize(n, 0);
            this->stampEarlier.resize(n, 0);
            this->stampIn.resize(n, 0);
            this->stampOut.resize(n, 0);
            this->stampNew.resize(n, 0);
        }
#ifndef BATCH_MINI_SCAN
        if ((int)this->minIbst.size() < n) {
            // empty trees ordered by the positions of this algorithm
            this->minIbst.resize(n, DifferentBalancedBinaryTree(this->reverseTopologicalOrder));
        }
#endif
    }


    // =======================================================================
    // Base graph
    // =======================================================================

    // the edge only goes into the graph; dfs() computes the order
    void IncrementalDFSTREEBatch::insert(GraphContext *g, int u, int v) {
        g->insertEdge(u, v);
        this->n = g->displayNumberVertices();
        this->fitSize(this->n);
    }


    // static DFS of the base graph (BaseDFS::dfs, as in static), then
    // IncrementalTopoOrder starts from exactly this order, with all edges of
    // the graph (including the artificial ones)
    void IncrementalDFSTREEBatch::dfs(GraphContext *g) {
        this->n = g->displayNumberVertices();
        this->fitSize(this->n);
        std::fill(this->visit.begin(), this->visit.end(), 0);
        this->postorderNumber = 0;
        BaseDFS::dfs(g);

        std::vector<inctopo::Edge> edges;
        for (int x = 0; x < this->n; x++) {
            for (int y : g->displayOutgoingNeighbours(x)) {
                edges.push_back(inctopo::Edge(x, y));
            }
        }
        std::vector<inctopo::Vertex> order(this->topOrd.begin(), this->topOrd.begin() + this->n);
        this->topo.assign(this->n, edges, order);
    }


    void IncrementalDFSTREEBatch::calculateMinI(GraphContext *g) {
        this->minI[0] = INT_MAX;
        this->minI[1] = 0;
#ifndef BATCH_MINI_SCAN
        // trees of the base graph: every vertex gets its real parents
        for (int i = 2; i < this->n; i++) {
            this->minIbst[i] = DifferentBalancedBinaryTree(this->reverseTopologicalOrder);
            for (int u : g->displayIncomingNeighbours(i)) {
                if (u != 0) this->minIbst[i].addNode(u, 0);
            }
            this->minI[i] = (g->strategy->incoming[i] == 2) ? 0 : this->minIbst[i].findMin();
        }
#else
        for (int i = 2; i < this->n; i++) this->minIOf(g, i);
#endif
    }


    void IncrementalDFSTREEBatch::calculateMaxO(GraphContext *g) {
        this->maxO[1] = INT_MAX;
        this->maxO[0] = 1;
        for (int i = 2; i < this->n; i++) this->maxOOf(g, i);
    }


    // =======================================================================
    // A batch
    // =======================================================================

    void IncrementalDFSTREEBatch::insertBatch(GraphContext *g) {
        Graph *G = g->strategy;
        this->beginBatch();
        BATCH_TIMER_START

        // ---- step 1: every edge of the batch
        for (const std::pair<int,int> &e : this->pending) {
            std::map<int,int>::iterator iu = G->reverseMapping.find(e.first);
            std::map<int,int>::iterator iv = G->reverseMapping.find(e.second);
            bool uExists = (iu != G->reverseMapping.end());
            bool vExists = (iv != G->reverseMapping.end());
            if (uExists && vExists) {
                const std::list<int> &out = G->adjLOutgoing[iu->second];
                if (std::find(out.begin(), out.end(), iv->second) != out.end()) {
                    this->duplicateEdges++;       // already in the graph: nothing changes
                    continue;
                }
            }
            // read before the graph changes its flags
            bool vWasSource = vExists && G->incoming[iv->second] == 2;
            int caseNo;
            if (!uExists && !vExists) caseNo = 3;
            else if (!uExists) caseNo = 1;
            else if (!vExists) caseNo = 2;
            else caseNo = (this->topo.position(iu->second) > this->topo.position(iv->second)) ? 4 : 5;

            g->insertEdge(e.first, e.second);
            this->fitSize(g->displayNumberVertices());
            int ru = G->reverseMapping[e.first], rv = G->reverseMapping[e.second];
            this->realEdges.push_back(std::make_pair(ru, rv));

            switch (caseNo) {
            case 1: this->edgeTypes.currentcase1++; this->placeNewVertex(1, ru, rv, vWasSource); break;
            case 2: this->edgeTypes.currentcase2++; this->placeNewVertex(2, ru, rv, false); break;
            case 3: this->edgeTypes.currentcase3++; this->placeNewVertex(3, ru, rv, false); break;
            case 4:
                // points backwards: repaired by the next reorder
                this->edgeTypes.currentcase4++;
                this->newEdges.push_back(inctopo::Edge(ru, rv));
                break;
            default:
                this->edgeTypes.currentcase5++;
                if (vWasSource) {
                    // adoption: the block must be read from a DFS order, so the
                    // pending reorder is done first; it may also have put u
                    // behind v, then the edge is an ordinary backward edge
                    this->flushReorder();
                    if (this->topo.position(ru) < this->topo.position(rv)) this->adoptSubtree(g, ru, rv);
                    else this->newEdges.push_back(inctopo::Edge(ru, rv));
                }
                else {
                    this->newEdges.push_back(inctopo::Edge(ru, rv));   // consistent; only into the adjacency
                }
            }
        }
        this->pending.clear();
        this->n = g->displayNumberVertices();
        this->vertexCount.currentTotalEdges = g->displayNumberEdges();

        // ---- step 2: one reorder for the rest of the batch
        this->flushReorder();
        BATCH_TIMER_ADD(this->latestTimes.currentinsertTime)
        // the reorders and block moves are timed separately (dfsTime)
        this->latestTimes.currentinsertTime -= this->latestTimes.currentdfsTime;

        // ---- step 3 and 4
#ifndef BATCH_MINI_SCAN
        this->detachMovedFromTrees(g);
#endif
        this->copyOrder();
        this->updateMinIMaxO(g);

        // ---- totals
        this->latestTimes.insertTime += this->latestTimes.currentinsertTime;
        this->latestTimes.dfsTime += this->latestTimes.currentdfsTime;
        this->latestTimes.topologicalOrderTime += this->latestTimes.currenttopologicalOrderTime;
        this->latestTimes.minITime += this->latestTimes.currentminITime;
        this->latestTimes.maxOTime += this->latestTimes.currentmaxOTime;
        this->edgeTypes.case1 += this->edgeTypes.currentcase1;
        this->edgeTypes.case2 += this->edgeTypes.currentcase2;
        this->edgeTypes.case3 += this->edgeTypes.currentcase3;
        this->edgeTypes.case4 += this->edgeTypes.currentcase4;
        this->edgeTypes.case5 += this->edgeTypes.currentcase5;
        this->vertexCount.countAffectedVertices += this->vertexCount.currentcountAffectedVertices;
        this->vertexCount.countExploredVertices += this->vertexCount.currentcountExploredVertices;
        this->vertexCount.totalEdges = this->vertexCount.currentTotalEdges;
    }


    // per-batch state: counters, lists, stamps (a new stamp value empties all
    // stamp sets at once)
    void IncrementalDFSTREEBatch::beginBatch() {
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

        if (++this->stamp == 0) {     // wrap-around: clear the stamps once
            std::fill(this->stampLater.begin(), this->stampLater.end(), 0);
            std::fill(this->stampEarlier.begin(), this->stampEarlier.end(), 0);
            std::fill(this->stampIn.begin(), this->stampIn.end(), 0);
            std::fill(this->stampOut.begin(), this->stampOut.end(), 0);
            this->stamp = 1;
        }
        this->newEdges.clear();
        this->backEdges.clear();
        this->realEdges.clear();
        this->movedLater.clear();
        this->movedEarlier.clear();
        this->dirtyLo = INT_MAX;
        this->dirtyHi = -1;
    }


    void IncrementalDFSTREEBatch::markDirty(int lo, int hi) {
        this->dirtyLo = std::min(this->dirtyLo, lo);
        this->dirtyHi = std::max(this->dirtyHi, hi);
    }


    /*
      Step 1, cases 1-3: a new vertex goes where the DFS would visit it (rules
      of incMinIBST); inserting it keeps the relative order of all others.
        case 1  u directly after 0, or directly before v if v's only parent
                was 0 (u takes v's place, v's subtree hangs under u)
        case 2  v directly after u
        case 3  u, v directly before 1
      Its edges (the real one and the artificial 0 -> u / v -> 1 of the graph)
      only go into the adjacency (none points backwards).
    */
    void IncrementalDFSTREEBatch::placeNewVertex(int caseNo, int ru, int rv, bool vWasSource) {
        int pos, a, b;
        bool ok = true;
        switch (caseNo) {
        case 1:
            pos = vWasSource ? this->topo.position(rv) : 1;
            a = this->topo.addVertexAt(pos);
            ok = (a == ru);
            this->newEdges.push_back(inctopo::Edge(0, ru));
            this->newEdges.push_back(inctopo::Edge(ru, rv));
            this->minI[ru] = 0;            // no real parent
            this->maxO[ru] = 1;            // set to v by updateMinIMaxO
            break;
        case 2:
            pos = this->topo.position(ru) + 1;
            b = this->topo.addVertexAt(pos);
            ok = (b == rv);
            this->newEdges.push_back(inctopo::Edge(ru, rv));
            this->newEdges.push_back(inctopo::Edge(rv, 1));
            this->minI[rv] = 0;            // set to u by updateMinIMaxO
            this->maxO[rv] = 1;            // no real child
            break;
        default:
            pos = this->topo.position(1);
            a = this->topo.addVertexAt(pos);
            b = this->topo.addVertexAt(this->topo.position(1));
            ok = (a == ru && b == rv);
            this->newEdges.push_back(inctopo::Edge(0, ru));
            this->newEdges.push_back(inctopo::Edge(ru, rv));
            this->newEdges.push_back(inctopo::Edge(rv, 1));
            this->minI[ru] = 0; this->maxO[ru] = 1;
            this->minI[rv] = 0; this->maxO[rv] = 1;
        }
        if (!ok) {
            std::cerr << "incDFSBatch: vertex ids of the graph and the order are out of step" << std::endl;
            std::exit(1);
        }
        // every position from pos on shifted by one
        this->markDirty(pos, this->topo.size() - 1);
        // new vertices are in no tree yet (and have no position in
        // reverseTopologicalOrder before step 3)
        if (caseNo != 2) this->stampNew[ru] = this->stamp;
        if (caseNo != 1) this->stampNew[rv] = this->stamp;
    }


    /*
      Step 1, adoption (case 5 with v a source; the edge u -> v is already in
      the graph, the order is a DFS order): v and its subtree must hang under
      u, as u's last child.
        - v's subtree is the block of positions from v up to the next child
          of 0 (a source, or the sink 1)
        - the block moves directly after u; the vertices in between move
          behind it
        - every edge from a vertex now behind the block into the block points
          backwards: that vertex is visited before u -> v in the DFS and
          takes over this part of v's old subtree; these edges are repaired by
          the next reorder, exactly like case 4 edges
    */
    void IncrementalDFSTREEBatch::adoptSubtree(GraphContext *g, int u, int v) {
        BATCH_TIMER_START
        Graph *G = g->strategy;
        const std::vector<int32_t> &inv = this->topo.order();
        int s = this->topo.position(v), e = s, last = this->topo.size() - 1;
        while (e + 1 <= last) {
            int x = inv[e + 1];
            if (x == 1 || (x >= 2 && G->incoming[x] == 2)) break;   // next child of 0
            e++;
        }
        int to = this->topo.position(u) + 1;
        this->newEdges.push_back(inctopo::Edge(u, v));      // into the adjacency
        if (to < s) {
            this->topo.moveBlock(s, e, to);
            this->markDirty(to, e);
            int blockEnd = to + (e - s);
            for (int p = to; p <= blockEnd; p++) {
                int w = inv[p];
                if (this->stampEarlier[w] != this->stamp) { this->stampEarlier[w] = this->stamp; this->movedEarlier.push_back(w); }
                for (int y : g->displayIncomingNeighbours(w)) {
                    if (y != 0 && this->topo.position(y) > p) {
                        this->backEdges.push_back(inctopo::Edge(y, w));
                        this->repairedEdges++;
                    }
                }
            }
        }
        this->adoptions++;
        BATCH_TIMER_ADD(this->latestTimes.currentdfsTime)
    }


    /*
      Step 2 (also before an adoption): one reorder for the pending new edges
      and the collected backward edges. The vertices it moves (reached by its
      searches) only move later relative to all others; they are recorded for
      minI / maxO.
    */
    void IncrementalDFSTREEBatch::flushReorder() {
        if (this->newEdges.empty() && this->backEdges.empty()) return;
        BATCH_TIMER_START
        inctopo::BatchStats st = this->topo.insertAndRepair(this->newEdges, this->backEdges);
        this->newEdges.clear();
        this->backEdges.clear();
        if (st.moved()) {
            this->markDirty(st.windowLo, st.windowHi);
            const std::vector<int32_t> &inv = this->topo.order();
            for (int p = st.windowLo; p <= st.windowHi; p++) {
                int x = inv[p];
                if (this->topo.movedInLastInsert(x) && this->stampLater[x] != this->stamp) {
                    this->stampLater[x] = this->stamp;
                    this->movedLater.push_back(x);
                }
            }
            this->vertexCount.currentcountExploredVertices += (int)st.dfsVisited;
        }
        BATCH_TIMER_ADD(this->latestTimes.currentdfsTime)
    }


    // Step 3: the changed positions into topOrd / reverseTopologicalOrder
    void IncrementalDFSTREEBatch::copyOrder() {
        BATCH_TIMER_START
        if (this->dirtyHi >= 0) {
            const std::vector<int32_t> &inv = this->topo.order();
            for (int p = this->dirtyLo; p <= this->dirtyHi; p++) {
                int x = inv[p];
                this->topOrd[p] = x;
                this->reverseTopologicalOrder[x] = p;
            }
            this->vertexCount.currentcountAffectedVertices = this->dirtyHi - this->dirtyLo + 1;
        }
        BATCH_TIMER_ADD(this->latestTimes.currenttopologicalOrderTime)
    }


    /*
      Step 4: minI / maxO (as in incMinIBST, extended for adoptions).
      The relative order of two vertices changes only if one of them
        - moved later in a reorder (movedLater), or
        - moved earlier in an adoption (movedEarlier, the moved blocks);
      inserting new vertices keeps every relative order. Hence:
        minI[c] (earliest parent)
          - rescan c's parents if its earliest parent moved later (another
            parent may now be the earliest) or if a parent moved earlier
          - compare with the batch's new parents of c
        maxO[q] (latest child)
          - rescan q's children if its latest child moved earlier
          - compare with q's children that moved later and its new children
            (a child that moved later can only become the latest)
      Old values are read before they are changed (the rescans come first;
      a rescanned vertex is final, it already sees all new edges).
    */
    void IncrementalDFSTREEBatch::updateMinIMaxO(GraphContext *g) {
#ifndef BATCH_MINI_SCAN
        {   // ---- minI (trees, as in incMinIBST)
        BATCH_TIMER_START
        this->attachMovedToTrees(g);
        BATCH_TIMER_ADD(this->latestTimes.currentminITime)
        }
#else
        {   // ---- minI
        BATCH_TIMER_START
        for (int x : this->movedLater) {
            for (int c : g->displayOutgoingNeighbours(x)) {
                if (c != 1 && this->minI[c] == x && this->stampIn[c] != this->stamp) {
                    this->stampIn[c] = this->stamp;
                    this->minIOf(g, c);
                }
            }
        }
        for (int w : this->movedEarlier) {
            for (int c : g->displayOutgoingNeighbours(w)) {
                if (c != 1 && this->stampIn[c] != this->stamp) {
                    this->stampIn[c] = this->stamp;
                    this->minIOf(g, c);
                }
            }
        }
        for (const std::pair<int,int> &e : this->realEdges) {
            int u = e.first, v = e.second;
            if (this->stampIn[v] == this->stamp) continue;
            if (this->minI[v] == 0 || this->reverseTopologicalOrder[u] < this->reverseTopologicalOrder[this->minI[v]]) {
                this->minI[v] = u;
            }
        }
        BATCH_TIMER_ADD(this->latestTimes.currentminITime)
        }
#endif
        {   // ---- maxO
        BATCH_TIMER_START
        for (int w : this->movedEarlier) {
            for (int q : g->displayIncomingNeighbours(w)) {
                if (q != 0 && this->maxO[q] == w && this->stampOut[q] != this->stamp) {
                    this->stampOut[q] = this->stamp;
                    this->maxOOf(g, q);
                }
            }
        }
        for (int x : this->movedLater) {
            for (int q : g->displayIncomingNeighbours(x)) {
                if (q == 0 || this->stampOut[q] == this->stamp) continue;
                if (this->maxO[q] == 1 || this->reverseTopologicalOrder[x] > this->reverseTopologicalOrder[this->maxO[q]]) {
                    this->maxO[q] = x;
                }
            }
        }
        for (const std::pair<int,int> &e : this->realEdges) {
            int u = e.first, v = e.second;
            if (this->stampOut[u] == this->stamp) continue;
            if (this->maxO[u] == 1 || this->reverseTopologicalOrder[v] > this->reverseTopologicalOrder[this->maxO[u]]) {
                this->maxO[u] = v;
            }
        }
        BATCH_TIMER_ADD(this->latestTimes.currentmaxOTime)
        }
    }


#ifndef BATCH_MINI_SCAN
    /*
      minI trees (as in incMinIBST's deleteReachableFromBST /
      insertReachableIntoBST). A tree compares its vertices by
      reverseTopologicalOrder; the vertices whose relative order changes are
      the moved ones (movedLater, movedEarlier), so they are removed from
      their children's trees while reverseTopologicalOrder still has the old
      positions (before step 3), and put back afterwards. New vertices are in
      no tree yet. Then the batch's new edges go into the trees, and minI is
      the minimum of every tree that changed.
    */
    void IncrementalDFSTREEBatch::detachMovedFromTrees(GraphContext *g) {
        BATCH_TIMER_START
        const std::vector<int> *lists[2] = { &this->movedLater, &this->movedEarlier };
        for (const std::vector<int> *list : lists) {
            for (int x : *list) {
                if (this->stampNew[x] == this->stamp) continue;
                for (int c : g->displayOutgoingNeighbours(x)) {
                    if (c != 1) this->minIbst[c].deleteNode(x);
                }
            }
        }
        BATCH_TIMER_ADD(this->latestTimes.currentminITime)
    }


    void IncrementalDFSTREEBatch::attachMovedToTrees(GraphContext *g) {
        std::vector<int> changed;
        const std::vector<int> *lists[2] = { &this->movedLater, &this->movedEarlier };
        for (const std::vector<int> *list : lists) {
            for (int x : *list) {
                for (int c : g->displayOutgoingNeighbours(x)) {
                    if (c == 1) continue;
                    this->minIbst[c].addNode(x, 0);
                    if (this->stampIn[c] != this->stamp) { this->stampIn[c] = this->stamp; changed.push_back(c); }
                }
            }
        }
        for (const std::pair<int,int> &e : this->realEdges) {
            this->minIbst[e.second].addNode(e.first, 0);
            if (this->stampIn[e.second] != this->stamp) { this->stampIn[e.second] = this->stamp; changed.push_back(e.second); }
        }
        for (int c : changed) {
            this->minI[c] = (g->strategy->incoming[c] == 2) ? 0 : this->minIbst[c].findMin();
        }
    }
#endif


    // =======================================================================
    // minI / maxO of one vertex, printing
    // =======================================================================

    // rules of SimpleDFS::calculateMinI: 0 if the only parent is the
    // artificial source, otherwise the earliest real parent
    void IncrementalDFSTREEBatch::minIOf(GraphContext *g, int x) {
        if (x < 2) return;
        if (g->strategy->incoming[x] == 2) {
            this->minI[x] = 0;
            return;
        }
        int best = INT_MAX, node = INT_MAX;
        for (int u : g->displayIncomingNeighbours(x)) {
            if (u == 0) continue;
            if (this->reverseTopologicalOrder[u] < best) {
                best = this->reverseTopologicalOrder[u];
                node = u;
            }
        }
        this->minI[x] = node;
    }


    // rules of SimpleDFS::calculateMaxO: 1 if the only child is the
    // artificial sink, otherwise the latest real child
    void IncrementalDFSTREEBatch::maxOOf(GraphContext *g, int x) {
        if (x < 2) return;
        if (g->strategy->outgoing[x] == 2) {
            this->maxO[x] = 1;
            return;
        }
        int best = -1, node = INT_MAX;
        for (int u : g->displayOutgoingNeighbours(x)) {
            if (u == 1) continue;
            if (this->reverseTopologicalOrder[u] > best) {
                best = this->reverseTopologicalOrder[u];
                node = u;
            }
        }
        this->maxO[x] = node;
    }


    void IncrementalDFSTREEBatch::printMinI() const {
        std::cout << "MinI of the graph is" << std::endl;
        for (int i = 0; i < this->n; i++) std::cout << i << " -> " << this->minI[i] << std::endl;
    }


    void IncrementalDFSTREEBatch::printMaxO() const {
        std::cout << "MaxO of the graph is " << std::endl;
        for (int i = 0; i < this->n; i++) std::cout << i << " -> " << this->maxO[i] << std::endl;
    }
