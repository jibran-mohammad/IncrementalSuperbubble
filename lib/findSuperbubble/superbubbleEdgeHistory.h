#ifndef SUPERBUBBLE_EDGE_HISTORY_H
#define SUPERBUBBLE_EDGE_HISTORY_H
#include <algorithm>
#include <fstream>
#include <iterator>
#include <cstdint>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>
#include "../graph/contextGraph.h"
#include "superbubbleTree.h"

/*
  ===========================================================================
  Superbubble history per edge: what every inserted edge changes
  ===========================================================================
  The per-edge counterpart of superbubbleHistory.h (which records per genome),
  for the heuristic incMinIBSTH1 in perEdge mode (third argument "history").
  The superbubble routine (SuperbubbleTree) journals every superbubble it
  creates and removes while processing one edge; only the NET change is
  recorded (a superbubble removed and created again by a rescan or a full
  recomputation of the same edge does not appear). A superbubble is
  identified by its (entrance, exit) pair; the artificial <0 1> is left out.

  Only superbubbles with at least one inside vertex are recorded (as in
  superbubbleHistory.h): a trivial superbubble (a single edge, nothing
  inside) is not a variant site. A pair that gains its first inside vertex
  counts as created, one that loses its last as destroyed. Every pair the
  routine touched is checked against the order after the edge (a
  superbubble is the range entrance .. exit of the DFS order, so it has an
  inside vertex iff the range is longer than the two).

  <prefix>_edge_events.tsv, one row per inserted edge (vertex ids are the
  ids of the input edge list):
    edge_idx       running number over all inserted edges (from 1)
    step           genome number (1 = first genome after the base graph)
    genome         name of the genome
    edge_in_genome number of the edge within its genome (from 1)
    u, v           the edge
    case           case of the incremental DFS (1-5); 0 = the edge was
                   already in the graph (nothing changes)
    full           1 if the routine recomputed all superbubbles
    n_destroyed    superbubbles (with an inside vertex) that are no longer
                   superbubbles or have nothing inside any more
    n_created      new superbubbles with an inside vertex (also a pair that
                   gains its first inside vertex)
    alive          superbubbles with an inside vertex after the edge
    destroyed, created
                   the superbubbles as entrance-exit, separated by ','
                   ("-": none)
    time_us        time of the edge (order, minI, maxO, superbubbles), us

  The superbubbles alive at the end are those of the base graph plus every
  created and minus every destroyed one, in edge order; they must equal the
  superbubbles with an inside vertex of the final superbubble list (a check
  on real data, Preprocessing/edgeHistory.py --final --graph).
  ===========================================================================
*/
class SuperbubbleEdgeHistory {
public:
    bool open(const std::string &prefix) {
        events.open(prefix + "_edge_events.tsv");
        if (!events.is_open()) return false;
        events << "edge_idx\tstep\tgenome\tedge_in_genome\tu\tv\tcase\tfull\tn_destroyed\tn_created\talive"
                  "\tdestroyed\tcreated\ttime_us\n";
        return true;
    }

    // after the base graph: the superbubbles with an inside vertex (pos: the
    // position of every vertex in the DFS order)
    void start(const SuperbubbleTree &tree, const std::vector<int> &pos) {
        nontrivial.clear();
        for (const std::pair<int,int> &sb : tree.superbubbles())
            if (pos[sb.second] - pos[sb.first] > 1) nontrivial.insert(key(sb.first, sb.second));
    }

    // after one edge: x, y the edge in input ids, caseNo 0 for a duplicate
    // edge; takes and clears the journal of the tree
    void record(int step, const std::string &genome, int edgeInGenome, int x, int y, int caseNo, bool full,
                SuperbubbleTree &tree, GraphContext &g, const std::vector<int> &pos, double timeSec) {
        std::vector<std::pair<int,int>> &touched = tree.journalAdded, &killed = tree.journalKilled;
        touched.insert(touched.end(), killed.begin(), killed.end());
        std::sort(touched.begin(), touched.end());
        touched.erase(std::unique(touched.begin(), touched.end()), touched.end());
        created.clear(); destroyed.clear();
        for (const std::pair<int,int> &sb : touched) {
            uint64_t k = key(sb.first, sb.second);
            bool before = nontrivial.count(k) > 0;
            bool after = tree.alive(sb.first, sb.second) && pos[sb.second] - pos[sb.first] > 1;
            if (before && !after) { destroyed.push_back(sb); nontrivial.erase(k); }
            if (!before && after) { created.push_back(sb); nontrivial.insert(k); }
        }
        touched.clear(); killed.clear();

        events << ++edgeIdx << '\t' << step << '\t' << genome << '\t' << edgeInGenome << '\t' << x << '\t' << y
               << '\t' << caseNo << '\t' << (full ? 1 : 0) << '\t' << destroyed.size() << '\t' << created.size()
               << '\t' << nontrivial.size() << '\t';
        list(destroyed, g);
        events << '\t';
        list(created, g);
        events << '\t' << (long long)(timeSec * 1e6 + 0.5) << '\n';
    }

private:
    std::ofstream events;
    long long edgeIdx = 0;
    std::vector<std::pair<int,int>> created, destroyed;
    std::unordered_set<uint64_t> nontrivial;          // current superbubbles with an inside vertex

    static uint64_t key(int s, int t) { return ((uint64_t)(uint32_t)s << 32) | (uint32_t)t; }

    void list(const std::vector<std::pair<int,int>> &sbs, GraphContext &g) {
        if (sbs.empty()) { events << '-'; return; }
        for (size_t i = 0; i < sbs.size(); ++i)
            events << (i ? "," : "") << g.toOriginal(sbs[i].first) << '-' << g.toOriginal(sbs[i].second);
    }
};
#endif
