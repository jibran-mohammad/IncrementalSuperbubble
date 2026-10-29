#include"superbubbleBatchH1.h"
#include<algorithm>
#include<chrono>
#include<iostream>
#include<numeric>
#include<unordered_map>
#include<unordered_set>

    // read-only view of the current order, minI and maxO (no copy)
    SBView SuperbubbleCalculationBatchH1::view(DFSCONTEXT *sdfs, GraphContext *g) {
        IDFS *d = sdfs->dfsStrategy;
        SBView V;
        V.n = g->displayNumberVertices();
        V.pos = d->reverseTopologicalOrder.data();
        V.at = d->topOrd.data();
        V.minI = d->minI.data();
        V.maxO = d->maxO.data();
        return V;
    }


    void SuperbubbleCalculationBatchH1::computeSB(DFSCONTEXT *sdfs, GraphContext *g) {
        this->tree.build(view(sdfs, g));
        this->pending = false;
    }


    void SuperbubbleCalculationBatchH1::insertBatch(DFSCONTEXT *sdfs, GraphContext *g, const std::vector<std::pair<int,int>> &edges) {
        #ifndef BENCHMARK
        auto start = std::chrono::high_resolution_clock::now();
        #endif
        this->readGroups(g, edges);
        #ifndef BENCHMARK
        this->prepareTime = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
        #endif

        sdfs->updateEdgeBatch(g, edges);
        this->pending = true;
    }


    /*
      The groups of a batch, from the graph before it. The graph marks a
      vertex whose only child is the sink 1 with outgoing == 2, and one whose
      only parent is the source 0 with incoming == 2.
      New vertices get local numbers; a union-find over the edges between two
      new vertices gives the groups. Per new vertex the in- and out-degree in
      the batch decide whether a dangling group is one simple path.
    */
    void SuperbubbleCalculationBatchH1::readGroups(GraphContext *g, const std::vector<std::pair<int,int>> &edges) {
        Graph *G = g->strategy;
        this->groups.clear();

        std::unordered_map<int,int> local;          // input id of a new vertex -> local number
        std::vector<int> input;                     // local number -> input id
        auto idOf = [&](int x) -> int {
            std::map<int,int>::const_iterator it = G->reverseMapping.find(x);
            return it == G->reverseMapping.end() ? -1 : it->second;
        };
        auto localOf = [&](int x) -> int {
            std::unordered_map<int,int>::const_iterator it = local.find(x);
            if (it != local.end()) return it->second;
            int k = (int)input.size();
            local[x] = k;
            input.push_back(x);
            return k;
        };

        // 1. edges: between two old vertices -> a group of its own; others kept
        struct NewEdge { int u, v; bool uNew, vNew; };   // old: internal id, new: local number
        std::vector<NewEdge> newEdges;
        std::unordered_set<unsigned long long> seen;
        for (const std::pair<int,int> &e : edges) {
            unsigned long long key = ((unsigned long long)(unsigned)e.first << 32) | (unsigned)e.second;
            if (!seen.insert(key).second) continue;           // twice in the batch
            int ru = idOf(e.first), rv = idOf(e.second);
            if (ru >= 0 && rv >= 0) {
                const std::list<int> &out = G->adjLOutgoing[ru];
                if (std::find(out.begin(), out.end(), rv) != out.end()) { this->duplicateEdges++; continue; }
                SuperbubbleTree::BatchGroup gr;
                gr.tails.push_back(ru);
                gr.heads.push_back(rv);
                gr.tailWasSink = (G->outgoing[ru] == 2);
                gr.headWasSource = (G->incoming[rv] == 2);
                this->groups.push_back(gr);
                continue;
            }
            NewEdge ne;
            ne.uNew = (ru < 0);
            ne.vNew = (rv < 0);
            ne.u = ne.uNew ? localOf(e.first) : ru;
            ne.v = ne.vNew ? localOf(e.second) : rv;
            newEdges.push_back(ne);
        }

        // 2. groups of new vertices (union-find)
        int m = (int)input.size();
        std::vector<int> parent(m);
        std::iota(parent.begin(), parent.end(), 0);
        auto find = [&](int x) -> int {
            while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
            return x;
        };
        for (const NewEdge &ne : newEdges) {
            if (ne.uNew && ne.vNew) parent[find(ne.u)] = find(ne.v);
        }
        int first = (int)this->groups.size();          // the groups of new vertices follow
        std::vector<int> groupOf(m), rootGroup(m, -1);
        std::vector<std::vector<int>> members;
        for (int k = 0; k < m; k++) {
            int r = find(k);
            if (rootGroup[r] < 0) { rootGroup[r] = (int)members.size(); members.push_back(std::vector<int>()); }
            groupOf[k] = rootGroup[r];
            members[groupOf[k]].push_back(k);
        }
        this->groups.resize(first + members.size());

        // 3. tails, heads, degrees, path successors
        std::vector<int> inDeg(m, 0), outDeg(m, 0), next(m, -1);
        std::vector<char> newParent(m, 0);
        std::unordered_map<int,int> newChildren, newParents;   // old vertex -> its edges to / from new vertices
        for (const NewEdge &ne : newEdges) {
            if (ne.uNew) outDeg[ne.u]++;
            if (ne.vNew) inDeg[ne.v]++;
            if (ne.uNew && ne.vNew) { next[ne.u] = ne.v; newParent[ne.v] = 1; }
            else if (ne.vNew) { this->groups[first + groupOf[ne.v]].tails.push_back(ne.u); newChildren[ne.u]++; }
            else { this->groups[first + groupOf[ne.u]].heads.push_back(ne.v); newParents[ne.v]++; }
        }

        for (int k = 0; k < (int)members.size(); k++) {
            SuperbubbleTree::BatchGroup &gr = this->groups[first + k];
            std::sort(gr.tails.begin(), gr.tails.end());
            gr.tails.erase(std::unique(gr.tails.begin(), gr.tails.end()), gr.tails.end());
            std::sort(gr.heads.begin(), gr.heads.end());
            gr.heads.erase(std::unique(gr.heads.begin(), gr.heads.end()), gr.heads.end());
            for (int t : gr.tails) if (G->outgoing[t] == 2) gr.tailWasSink = true;
            for (int h : gr.heads) if (G->incoming[h] == 2) gr.headWasSource = true;
            if (!gr.tails.empty() && !gr.heads.empty()) continue;

            // dangling: one simple path (at most one tail and one head, every
            // new vertex at most one parent and one child, one start)
            bool path = gr.tails.size() <= 1 && gr.heads.size() <= 1;
            int startV = -1;
            for (int x : members[k]) {
                if (inDeg[x] > 1 || outDeg[x] > 1) path = false;
                if (!newParent[x]) { if (startV >= 0) path = false; startV = x; }
            }
            if (!path || startV < 0) continue;
            for (int x = startV; x >= 0; x = next[x]) gr.path.push_back(input[x]);
            if (gr.path.size() != members[k].size()) { gr.path.clear(); continue; }
            gr.isPath = true;
            gr.tailPair = !gr.tails.empty() && gr.tailWasSink && newChildren[gr.tails[0]] == 1;
            gr.headPair = !gr.heads.empty() && gr.headWasSource && newParents[gr.heads[0]] == 1;
        }
    }


    void SuperbubbleCalculationBatchH1::updateSB(DFSCONTEXT *sdfs, GraphContext *g) {
        if (!this->pending) return;
        this->pending = false;
        // the new vertices of the paths: input ids -> internal ids
        Graph *G = g->strategy;
        for (SuperbubbleTree::BatchGroup &gr : this->groups) {
            for (int &x : gr.path) x = G->reverseMapping.find(x)->second;
        }
        this->tree.updateBatch(view(sdfs, g), this->groups);
    }


    std::vector<std::pair<int,int>> SuperbubbleCalculationBatchH1::superbubbles() const {
        std::vector<std::pair<int,int>> out = this->tree.superbubbles();
        out.push_back(std::make_pair(0, 1));
        return out;
    }


    // internal vertex ids are converted to input ids only for printing
    void SuperbubbleCalculationBatchH1::printSuperbubbles(GraphContext *g) {
    	std::cout << "superbubbles of graph are" << std::endl;
    	for (const std::pair<int,int> &sb : this->superbubbles()) {
            std::cout << "<" << g->toOriginal(sb.first) << " " << g->toOriginal(sb.second) << ">" << std::endl;
        }
        std::cout << std::endl;
    }


    // One superbubble per line, "<entrance> <exit>", in input ids.
    void SuperbubbleCalculationBatchH1::writeSuperbubbles(GraphContext *g, std::ostream &out) {
        for (const std::pair<int,int> &sb : this->superbubbles()) {
            out << g->toOriginal(sb.first) << " " << g->toOriginal(sb.second) << "\n";
        }
    }


    int SuperbubbleCalculationBatchH1::numberOfSuperbubbles() {
    	return this->tree.count() + 1;   // + <0 1>
    }
