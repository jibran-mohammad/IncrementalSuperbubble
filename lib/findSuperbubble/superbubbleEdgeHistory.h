#ifndef SUPERBUBBLE_EDGE_HISTORY_H
#define SUPERBUBBLE_EDGE_HISTORY_H
#include <algorithm>
#include <fstream>
#include <iterator>
#include <string>
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
    n_destroyed    superbubbles that are no longer superbubbles
    n_created      new superbubbles
    alive          superbubbles after the edge
    destroyed, created
                   the superbubbles as entrance-exit, separated by ','
                   ("-": none)
    time_us        time of the edge (order, minI, maxO, superbubbles), us

  The superbubbles alive at the end are those of the base graph plus every
  created and minus every destroyed one, in edge order; they must equal the
  final superbubble list (a check on real data).
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

    // after one edge: x, y the edge in input ids, caseNo 0 for a duplicate
    // edge; takes and clears the journal of the tree
    void record(int step, const std::string &genome, int edgeInGenome, int x, int y, int caseNo, bool full,
                SuperbubbleTree &tree, GraphContext &g, double timeSec) {
        std::vector<std::pair<int,int>> &added = tree.journalAdded, &killed = tree.journalKilled;
        std::sort(added.begin(), added.end());
        std::sort(killed.begin(), killed.end());
        created.clear(); destroyed.clear();
        std::set_difference(added.begin(), added.end(), killed.begin(), killed.end(), std::back_inserter(created));
        std::set_difference(killed.begin(), killed.end(), added.begin(), added.end(), std::back_inserter(destroyed));
        added.clear(); killed.clear();

        events << ++edgeIdx << '\t' << step << '\t' << genome << '\t' << edgeInGenome << '\t' << x << '\t' << y
               << '\t' << caseNo << '\t' << (full ? 1 : 0) << '\t' << destroyed.size() << '\t' << created.size()
               << '\t' << tree.count() << '\t';
        list(destroyed, g);
        events << '\t';
        list(created, g);
        events << '\t' << (long long)(timeSec * 1e6 + 0.5) << '\n';
    }

private:
    std::ofstream events;
    long long edgeIdx = 0;
    std::vector<std::pair<int,int>> created, destroyed;

    void list(const std::vector<std::pair<int,int>> &sbs, GraphContext &g) {
        if (sbs.empty()) { events << '-'; return; }
        for (size_t i = 0; i < sbs.size(); ++i)
            events << (i ? "," : "") << g.toOriginal(sbs[i].first) << '-' << g.toOriginal(sbs[i].second);
    }
};
#endif
