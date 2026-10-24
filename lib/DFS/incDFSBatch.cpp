#include"incDFSBatch.h"
#include<algorithm>
#include<cstdlib>
#include<unordered_map>
#include<unordered_set>

/*
  Batch algorithm: see the description in incDFSBatch.h. Order of this file:
    base graph          insert, dfs, calculateMinI / calculateMaxO
    a batch             insertBatch (steps 1-4) and its helpers in step order:
                        beginBatch, noteVertex / materialize (steps 1-2),
                        repairRound / search, adoptSubtree (step 3),
                        updateMinIMaxO (step 4)
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
            this->stampParent.resize(n, 0);
            this->stampDetached.resize(n, 0);
            this->markEpoch.resize(n, 0);
            this->tailEpoch.resize(n, 0);
            this->tailFirst.resize(n, -1);
            this->tailLast.resize(n, -1);
            this->afterList.resize(n);
            this->beforeList.resize(n);
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

    void IncrementalDFSTREEBatch::insert(GraphContext *g, int u, int v) {
        g->insertEdge(u, v);
        this->n = g->displayNumberVertices();
        this->fitSize(this->n);
    }


    // static DFS of the base graph (BaseDFS::dfs, as in the static algorithm):
    // topOrd / reverseTopologicalOrder
    void IncrementalDFSTREEBatch::dfs(GraphContext *g) {
        this->n = g->displayNumberVertices();
        this->fitSize(this->n);
        std::fill(this->visit.begin(), this->visit.end(), 0);
        this->postorderNumber = 0;
        BaseDFS::dfs(g);
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

    void IncrementalDFSTREEBatch::insertBatch(GraphContext *g, const std::vector<std::pair<int,int>> &edges) {
        Graph *G = g->strategy;
        this->beginBatch();
        const int nBefore = g->displayNumberVertices();
        std::unordered_map<int,int> newId;                  // input id of a new vertex -> its internal id
        int nAfter;

        {   // ---- step 1: classify (the graph is not changed)
        BATCH_TIMER_START
        std::unordered_set<unsigned long long> batchEdges;  // edges seen in this batch (input ids)
        auto known = [&](int x, int &id) -> bool {
            std::map<int,int>::const_iterator it = G->reverseMapping.find(x);
            if (it != G->reverseMapping.end()) { id = it->second; return true; }
            std::unordered_map<int,int>::const_iterator jt = newId.find(x);
            if (jt != newId.end()) { id = jt->second; return true; }
            return false;
        };
        // the graph numbers new vertices consecutively, in order of first appearance
        auto fresh = [&](int x) -> int {
            int id = nBefore + (int)newId.size();
            newId[x] = id;
            this->fitSize(id + 1);
            this->stampNew[id] = this->stamp;
            this->minI[id] = 0;              // no real parent yet
            this->maxO[id] = 1;              // no real child yet
            return id;
        };
        // a real parent in the graph, or one from an earlier edge of the batch
        auto hasParent = [&](int id) -> bool {
            return (id < nBefore && G->incoming[id] == 1) || this->stampParent[id] == this->stamp;
        };
        // A new vertex that is a source when it appears (u of case 1 or 3) is
        // placed only once the pass knows whether the genome gives it a parent:
        // then it goes directly after its first parent (where an adoption or a
        // reorder would put it), otherwise to its default place, noted at the
        // end in the order the sources appeared (the default lists hold only
        // such sources, so their order is that of noting them one by one).
        std::vector<int> dAnchor, dSeq;              // per new vertex (id - nBefore)
        std::vector<char> dBefore, dDeferred;
        int seq = 0;
        auto defer = [&](int s, int anchor, bool before, int sq) {
            int k = s - nBefore;
            if ((int)dDeferred.size() <= k) { dAnchor.resize(k + 1); dSeq.resize(k + 1); dBefore.resize(k + 1); dDeferred.resize(k + 1, 0); }
            dAnchor[k] = anchor; dBefore[k] = before; dSeq[k] = sq; dDeferred[k] = 1;
        };
        auto deferred = [&](int id) -> bool {
            int k = id - nBefore;
            return k >= 0 && k < (int)dDeferred.size() && dDeferred[k];
        };

        for (const std::pair<int,int> &e : edges) {
            int ru = -1, rv = -1;
            bool ku = known(e.first, ru), kv = known(e.second, rv);
            unsigned long long key = ((unsigned long long)(unsigned)e.first << 32) | (unsigned)e.second;
            if (ku && kv) {
                bool dup = batchEdges.count(key) > 0;
                if (!dup && ru < nBefore && rv < nBefore) {
                    const std::list<int> &adj = G->adjLOutgoing[ru];
                    dup = std::find(adj.begin(), adj.end(), rv) != adj.end();
                }
                if (dup) { this->duplicateEdges++; continue; }
                batchEdges.insert(key);
                // only a source whose subtree is made of this genome's new
                // vertices (default place after 0 or before the sink); a source
                // in an old source's place (case 1, v's only parent 0) carries
                // that old subtree, so its parent edge is an adoption (E_old)
                int k = rv - nBefore;
                if (deferred(rv) && !hasParent(rv) && (dAnchor[k] == 0 || dAnchor[k] == 1)) {
                    // the first parent of a source of this genome: v goes
                    // directly after u (instead of an adoption or a reorder)
                    dDeferred[rv - nBefore] = 0;
                    this->noteVertex(rv, ru, false);
                    this->sourcesPlaced++;
                    this->edgeTypes.currentcase5++;
                    this->newInputEdges.push_back(e);
                    this->realEdges.push_back(std::make_pair(ru, rv));
                }
                else {
                    this->oldEdges.push_back(OldEdge{e.first, e.second, ru, rv, !hasParent(rv)});
                }
            }
            else {
                batchEdges.insert(key);
                if (!ku && !kv) {
                    // case 3: u, v directly before the sink (u a source, v its child)
                    ru = fresh(e.first); rv = fresh(e.second);
                    this->edgeTypes.currentcase3++;
                    defer(ru, 1, true, seq++);
                    this->noteVertex(rv, ru, false);
                }
                else if (!ku) {
                    // case 1: u directly after 0, or in v's place if v's only parent was 0
                    ru = fresh(e.first);
                    this->edgeTypes.currentcase1++;
                    if (hasParent(rv)) defer(ru, 0, false, seq++);
                    else if (deferred(rv)) {
                        // v is a source of this genome: u takes its (default) place, v hangs under u
                        int k = rv - nBefore;
                        defer(ru, dAnchor[k], dBefore[k], dSeq[k]);
                        dDeferred[k] = 0;
                        this->noteVertex(rv, ru, false);
                    }
                    else defer(ru, rv, true, seq++);
                }
                else {
                    // case 2: v directly after u (u may itself be new)
                    rv = fresh(e.second);
                    this->edgeTypes.currentcase2++;
                    this->noteVertex(rv, ru, false);
                }
                this->newInputEdges.push_back(e);
                this->realEdges.push_back(std::make_pair(ru, rv));
            }
            this->stampParent[rv] = this->stamp;
        }
        nAfter = nBefore + (int)newId.size();
        {   // the sources without a parent in this genome: at their default places
            std::vector<int> rest;
            for (int k = 0; k < (int)dDeferred.size(); k++) if (dDeferred[k]) rest.push_back(nBefore + k);
            std::sort(rest.begin(), rest.end(), [&](int a, int b) { return dSeq[a - nBefore] < dSeq[b - nBefore]; });
            for (int s : rest) this->noteVertex(s, dAnchor[s - nBefore], dBefore[s - nBefore]);
        }

        // E_new into the graph, in order: the graph gives the new vertices
        // exactly the ids of step 1
        for (const std::pair<int,int> &e : this->newInputEdges) g->insertEdge(e.first, e.second);
        for (const std::pair<const int,int> &p : newId) {
            if (G->reverseMapping[p.first] != p.second) {
                std::cerr << "incDFSBatch: vertex ids of the graph and the order are out of step" << std::endl;
                std::exit(1);
            }
        }
        this->n = nAfter;
        this->fitSize(this->n);
        BATCH_TIMER_ADD(this->latestTimes.currentinsertTime)
        }

        // ---- step 2: build the new vertices into the order
        this->materialize(nBefore, nAfter);
        // a source placed under its parent may have an edge to an old vertex
        // before that parent: such E_new edges point backwards and go into the
        // first round (realEdges holds only E_new so far)
        for (const std::pair<int,int> &e : this->realEdges) {
            if (this->reverseTopologicalOrder[e.first] > this->reverseTopologicalOrder[e.second]) this->pending.push_back(e);
        }

        // ---- step 3: E_old, in order
        // During a round every backward edge of the graph must be one of the
        // round's edges (its searches follow all edges of the graph), so an
        // edge goes into the graph only when it is collected or consistent.
        for (const OldEdge &e : this->oldEdges) {
            // adoption (v's only parent was 0 and u is before v): the block
            // must be read from a DFS order, so the collected backward edges
            // are repaired first (before this edge is in the graph); they may
            // also have put u behind v, then the edge is an ordinary backward
            // edge. If u is already behind v, it is a backward edge anyway.
            bool adoption = e.vOnly0 && this->reverseTopologicalOrder[e.u] < this->reverseTopologicalOrder[e.v];
            if (adoption) this->repairRound(g);
            {
            BATCH_TIMER_START
            g->insertEdge(e.x, e.y);
            this->realEdges.push_back(std::make_pair(e.u, e.v));
            BATCH_TIMER_ADD(this->latestTimes.currentinsertTime)
            }
            if (adoption) {
                if (this->reverseTopologicalOrder[e.u] < this->reverseTopologicalOrder[e.v]) {
                    this->edgeTypes.currentcase5++;
                    this->adoptSubtree(g, e.u, e.v);
                }
                else {
                    this->edgeTypes.currentcase4++;
                    this->pending.push_back(std::make_pair(e.u, e.v));
                }
            }
            else if (this->reverseTopologicalOrder[e.u] > this->reverseTopologicalOrder[e.v]) {
                this->edgeTypes.currentcase4++;              // points backwards
                this->pending.push_back(std::make_pair(e.u, e.v));
            }
            else {
                this->edgeTypes.currentcase5++;              // the order stays valid
            }
        }
        this->repairRound(g);
        this->vertexCount.currentTotalEdges = g->displayNumberEdges();

        // ---- step 4: minI / maxO
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
        this->vertexCount.totalAffectedEdges += this->vertexCount.currentTotalAffectedEdges;
        this->vertexCount.totalExploredEdges += this->vertexCount.currentTotalExploredEdges;
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
            for (std::vector<unsigned> *s : { &this->stampLater, &this->stampEarlier, &this->stampIn, &this->stampOut,
                                              &this->stampNew, &this->stampParent, &this->stampDetached }) {
                std::fill(s->begin(), s->end(), 0);
            }
            this->stamp = 1;
        }
        this->newInputEdges.clear();
        this->oldEdges.clear();
        this->pending.clear();
        this->realEdges.clear();
        this->movedLater.clear();
        this->movedEarlier.clear();
        this->detached.clear();
    }


    // Step 1: a new vertex is noted next to its anchor (before or after it);
    // materialize() builds the noted vertices into the order in step 2
    void IncrementalDFSTREEBatch::noteVertex(int x, int anchor, bool before) {
        if (this->afterList[anchor].empty() && this->beforeList[anchor].empty()) this->anchors.push_back(anchor);
        (before ? this->beforeList[anchor] : this->afterList[anchor]).push_back(x);
    }


    /*
      Step 2: builds all new vertices into the order in one pass over it:
      every vertex x is emitted as  (its before-list, in noting order), x,
      (its after-list, newest first); a listed vertex is emitted the same way
      (chains of new vertices). Every new vertex is in exactly one list (that
      of its anchor), so every vertex is emitted once: O(n + new vertices).
      This is the order that inserting them one at a time gives. The new
      order is written into topOrd and reverseTopologicalOrder from the first
      new vertex on (the positions before it do not change).
    */
    void IncrementalDFSTREEBatch::materialize(int nBefore, int nAfter) {
        if (nAfter == nBefore) return;
        BATCH_TIMER_START
        std::vector<int> order;
        order.reserve(nAfter);
        int firstNew = nAfter;
        struct Frame { int v; int phase; int i; };
        std::vector<Frame> st;
        for (int p = 0; p < nBefore; p++) {
            st.push_back({this->topOrd[p], 0, 0});
            while (!st.empty()) {
                Frame &f = st.back();
                const std::vector<int> &before = this->beforeList[f.v];
                const std::vector<int> &after = this->afterList[f.v];
                if (f.phase == 0) {                                // before-list, in noting order
                    if (f.i < (int)before.size()) { int c = before[f.i++]; st.push_back({c, 0, 0}); continue; }
                    if (f.v >= nBefore && firstNew == nAfter) firstNew = (int)order.size();
                    order.push_back(f.v);                          // the vertex itself
                    f.phase = 1;
                    f.i = (int)after.size() - 1;
                }
                if (f.i >= 0) { int c = after[f.i--]; st.push_back({c, 0, 0}); continue; }   // after-list, newest first
                st.pop_back();
            }
        }
        if ((int)order.size() != nAfter) {
            std::cerr << "incDFSBatch: building in the new vertices lost some (" << order.size() << " of " << nAfter << ")" << std::endl;
            std::exit(1);
        }
        for (int p = firstNew; p < nAfter; p++) {
            this->topOrd[p] = order[p];
            this->reverseTopologicalOrder[order[p]] = p;
        }
        for (int a : this->anchors) { this->afterList[a].clear(); this->beforeList[a].clear(); }
        this->anchors.clear();
        BATCH_TIMER_ADD(this->latestTimes.currenttopologicalOrderTime)
    }


    /*
      Step 3: one reorder round for the collected backward edges (the batch
      DFS algorithm, IncrementalDFS + Update):
        window   lo = earliest head, hi = latest tail of the backward edges
        scan     positions hi .. lo (right to left); a vertex already reached
                 by a search is skipped (it has its place behind its tail);
                 for a tail x, for each of its backward edges (x, v) in their
                 order: search(v) collects, in postorder, the vertices reached
                 from v that lie before x (they must come behind x); then x
        relabel  the emitted vertices, read from right to left, are the new
                 window: one pass
      Each vertex of the window is emitted exactly once. The vertices that a
      search reached only move later relative to the others; they are taken
      out of the minI trees before the relabel and recorded for step 4.
    */
    void IncrementalDFSTREEBatch::repairRound(GraphContext *g) {
        if (this->pending.empty()) return;
        long long explored = 0;
        int lo = INT_MAX, hi = -1;
        {
        BATCH_TIMER_START
        if (++this->epoch == 0) {
            std::fill(this->markEpoch.begin(), this->markEpoch.end(), 0);
            std::fill(this->tailEpoch.begin(), this->tailEpoch.end(), 0);
            this->epoch = 1;
        }
        // the edges that still point backwards, per tail in their order
        this->edgeNext.assign(this->pending.size(), -1);
        for (int i = 0; i < (int)this->pending.size(); i++) {
            int u = this->pending[i].first, v = this->pending[i].second;
            if (this->reverseTopologicalOrder[u] < this->reverseTopologicalOrder[v]) continue;   // already consistent
            if (this->tailEpoch[u] != this->epoch) { this->tailEpoch[u] = this->epoch; this->tailFirst[u] = i; }
            else this->edgeNext[this->tailLast[u]] = i;
            this->tailLast[u] = i;
            lo = std::min(lo, this->reverseTopologicalOrder[v]);
            hi = std::max(hi, this->reverseTopologicalOrder[u]);
        }
        if (hi >= 0) {
            this->rounds++;
            this->out.clear();
            for (int p = hi; p >= lo; p--) {
                int x = this->topOrd[p];
                if (this->markEpoch[x] == this->epoch) continue;            // already behind its tail
                if (this->tailEpoch[x] == this->epoch) {
                    int bound = this->reverseTopologicalOrder[x];
                    for (int i = this->tailFirst[x]; i != -1; i = this->edgeNext[i]) {
                        explored += this->search(g, this->pending[i].second, x, bound);
                    }
                }
                this->out.push_back(x);
            }
        }
        BATCH_TIMER_ADD(this->latestTimes.currentdfsTime)
        }
        this->pending.clear();
        if (hi < 0) return;

        // the searched vertices change their relative order: out of the
        // minI trees while the positions are still the old ones
        for (int x : this->out) {
            if (this->markEpoch[x] == this->epoch) this->detachFromTrees(g, x);
        }
        {
        BATCH_TIMER_START
        int pos = hi;
        for (int x : this->out) {
            this->topOrd[pos] = x;
            this->reverseTopologicalOrder[x] = pos;
            pos--;
        }
        this->vertexCount.currentcountAffectedVertices += hi - lo + 1;
        this->vertexCount.currentcountExploredVertices += (int)explored;
        BATCH_TIMER_ADD(this->latestTimes.currenttopologicalOrderTime)
        }
    }


    /*
      The search of a reorder round (IncrementalDFS of the batch DFS
      algorithm), iterative: from v over the outgoing edges, only to
      vertices before the tail (position < bound) not yet reached in this
      round; a vertex is emitted when all its edges are done (postorder), so
      the deepest vertex comes last in the window. Reaching the tail itself
      would mean a cycle.
    */
    long long IncrementalDFSTREEBatch::search(GraphContext *g, int v, int tail, int bound) {
        if (this->markEpoch[v] == this->epoch || this->reverseTopologicalOrder[v] >= bound) return 0;
        struct Frame { int v; std::list<int>::const_iterator it; };
        std::vector<Frame> st;
        long long cnt = 1;
        this->markEpoch[v] = this->epoch;
        st.push_back({v, g->displayOutgoingNeighbours(v).begin()});
        while (!st.empty()) {
            Frame &f = st.back();
            const std::list<int> &adj = g->displayOutgoingNeighbours(f.v);
            if (f.it != adj.end()) {
                int w = *f.it;
                ++f.it;
                if (w == tail) {
                    std::cerr << "incDFSBatch: the batch creates a cycle" << std::endl;
                    std::exit(1);
                }
                if (this->markEpoch[w] == this->epoch || this->reverseTopologicalOrder[w] >= bound) continue;
                this->markEpoch[w] = this->epoch;
                cnt++;
                st.push_back({w, g->displayOutgoingNeighbours(w).begin()});
            }
            else {
                int x = f.v;
                st.pop_back();
                this->out.push_back(x);
                this->vertexCount.currentTotalExploredEdges += g->getOutDegree(x);
                if (this->stampLater[x] != this->stamp) { this->stampLater[x] = this->stamp; this->movedLater.push_back(x); }
            }
        }
        return cnt;
    }


    /*
      Step 3, adoption (case 5, v's only parent was 0; the edge u -> v is in
      the graph, the order is a DFS order): v and its subtree must hang under
      u, as u's last child.
        - v's subtree is the block of positions from v up to the next child
          of 0 (a source, or the sink 1)
        - the block moves directly after u; the vertices in between move
          behind it
        - every edge from a vertex now behind the block into the block points
          backwards: that vertex is visited before u -> v in the DFS and
          takes over this part of v's old subtree; these edges are repaired
          by the next round, like case 4 edges
    */
    void IncrementalDFSTREEBatch::adoptSubtree(GraphContext *g, int u, int v) {
        Graph *G = g->strategy;
        int s = this->reverseTopologicalOrder[v], e = s, last = this->n - 1;
        while (e + 1 <= last) {
            int x = this->topOrd[e + 1];
            if (x == 1 || (x >= 2 && G->incoming[x] == 2)) break;   // next child of 0
            e++;
        }
        int to = this->reverseTopologicalOrder[u] + 1;
        this->adoptions++;
        if (to >= s) return;
        // the block changes its relative order: out of the minI trees first
        for (int p = s; p <= e; p++) this->detachFromTrees(g, this->topOrd[p]);
        BATCH_TIMER_START
        std::rotate(this->topOrd.begin() + to, this->topOrd.begin() + s, this->topOrd.begin() + e + 1);
        for (int p = to; p <= e; p++) this->reverseTopologicalOrder[this->topOrd[p]] = p;
        this->vertexCount.currentcountAffectedVertices += e - to + 1;
        int blockEnd = to + (e - s);
        for (int p = to; p <= blockEnd; p++) {
            int w = this->topOrd[p];
            if (this->stampEarlier[w] != this->stamp) { this->stampEarlier[w] = this->stamp; this->movedEarlier.push_back(w); }
            for (int y : g->displayIncomingNeighbours(w)) {
                if (y != 0 && this->reverseTopologicalOrder[y] > p) {
                    this->pending.push_back(std::make_pair(y, w));
                    this->repairedEdges++;
                }
            }
        }
        BATCH_TIMER_ADD(this->latestTimes.currentdfsTime)
    }


    /*
      minI trees (as in incMinIBST's deleteReachableFromBST /
      insertReachableIntoBST). A tree compares its vertices by
      reverseTopologicalOrder, so a vertex is taken out of its children's
      trees before its relative order changes (once per batch) and put back
      in step 4. New vertices are in no tree yet.
    */
    void IncrementalDFSTREEBatch::detachFromTrees(GraphContext *g, int x) {
#ifndef BATCH_MINI_SCAN
        if (this->stampNew[x] == this->stamp || this->stampDetached[x] == this->stamp) return;
        BATCH_TIMER_START
        this->stampDetached[x] = this->stamp;
        this->detached.push_back(x);
        for (int c : g->displayOutgoingNeighbours(x)) {
            if (c != 1) this->minIbst[c].deleteNode(x);
        }
        BATCH_TIMER_ADD(this->latestTimes.currentminITime)
#endif
    }


    /*
      Step 4: minI / maxO (as in incMinIBST, extended for adoptions).
      The relative order of two vertices changes only if one of them
        - moved later in a reorder round (movedLater), or
        - moved earlier in an adoption (movedEarlier, the moved blocks);
      building in new vertices keeps every relative order. Hence:
        minI[c] (earliest parent)
          - trees: the detached vertices and the batch's new edges go into
            the trees, minI of every changed tree is its minimum
          - without trees: rescan c's parents if its earliest parent moved
            later or a parent moved earlier; compare with c's new parents
        maxO[q] (latest child)
          - rescan q's children if its latest child moved earlier
          - compare with q's children that moved later and its new children
            (a child that moved later can only become the latest)
    */
    void IncrementalDFSTREEBatch::updateMinIMaxO(GraphContext *g) {
#ifndef BATCH_MINI_SCAN
        {   // ---- minI (trees, as in incMinIBST)
        BATCH_TIMER_START
        this->attachToTrees(g);
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
    // the detached vertices back into their children's trees (a vertex that
    // is already in a tree is not added twice), then the batch's new edges;
    // minI of every changed tree is its minimum
    void IncrementalDFSTREEBatch::attachToTrees(GraphContext *g) {
        std::vector<int> changed;
        for (int x : this->detached) {
            for (int c : g->displayOutgoingNeighbours(x)) {
                if (c == 1) continue;
                this->minIbst[c].addNode(x, 0);
                if (this->stampIn[c] != this->stamp) { this->stampIn[c] = this->stamp; changed.push_back(c); }
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
