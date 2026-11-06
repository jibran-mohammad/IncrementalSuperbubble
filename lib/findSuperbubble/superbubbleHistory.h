#ifndef SUPERBUBBLE_HISTORY_H
#define SUPERBUBBLE_HISTORY_H
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iterator>
#include <list>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include "../graph/contextGraph.h"
#include "../DFS/contextdfs.h"

/*
  ===========================================================================
  Superbubble history: what every genome (batch) changes
  ===========================================================================
  After the base graph and after every genome, the current superbubble set is
  compared with the previous one, so only the NET change of a genome is
  recorded (superbubbles created and destroyed again within one genome, an
  artifact of the edge order, do not appear). A superbubble is identified by
  its (entrance, exit) pair; the artificial <0 1> is left out.

  Only superbubbles with at least one inside vertex are recorded. A trivial
  superbubble (a single edge s -> t, nothing inside) is not a variant site:
  the base graph (the reference path) consists of them only, one per edge.
  A trivial pair that gains an inside vertex is recorded as born / promoted
  at that step, one that loses all of them as died. The trivial ones are
  only counted (summary column "trivial").

  Two files (tab separated, vertex ids are the ids of the input edge list):

  <prefix>_history_events.tsv, one row per change:
    step       0 = base graph, k = after the k-th genome
    genome     name of the genome (base graph: "base")
    event      base     a superbubble of the base graph
               born     a new superbubble containing a vertex this genome
                        added (a new site)
               promoted a new superbubble made only of vertices that existed
                        before (an existing region became a superbubble)
               died     a superbubble that is no longer one
    entrance, exit
    n_inner    vertices strictly inside (base/born/promoted: now; died: when
               it was created)
    parent     base/born/promoted: entrance-exit of the innermost current
               superbubble containing it ("-": top level);
               died: entrance-exit of the innermost current superbubble that
               contains its former entrance and exit (absorbed into it; "-":
               none, top level)

  <prefix>_history_summary.tsv, one row per step:
    step genome new_edges new_vertices born promoted died died_merged alive time_s trivial
    alive = superbubbles with an inside vertex; trivial = the others (alive +
    trivial = all superbubbles, as in the time file without <0 1>)
    died_merged = died superbubbles absorbed by one born or promoted in the
    same step (smaller sites merged into a larger one)

  The superbubbles alive after step k are all rows created at a step <= k and
  not died at a step <= k; after the last step they must equal the final
  superbubble list (a check on real data).

  Membership uses the order: in a DFS order a superbubble <s t> is the
  interval of positions rto[s] .. rto[t]. New vertices of a step have the
  internal ids >= the vertex count before the step.
  ===========================================================================
*/
class SuperbubbleHistory {
public:
    bool open(const std::string &prefix) {
        events.open(prefix + "_history_events.tsv");
        summary.open(prefix + "_history_summary.tsv");
        if (!events.is_open() || !summary.is_open()) return false;
        events << "step\tgenome\tevent\tentrance\texit\tn_inner\tparent\n";
        summary << "step\tgenome\tnew_edges\tnew_vertices\tborn\tpromoted\tdied\tdied_merged\talive\ttime_s\ttrivial\n";
        return true;
    }

    // after the base graph (step 0) or a genome (step >= 1): SB is the current
    // superbubble list (internal ids, may contain <0 1>), nBefore the vertex
    // count before this step (0 for the base graph)
    void record(int step, const std::string &genome, DFSCONTEXT &d, GraphContext &g,
                const std::list<std::pair<int,int>> &SB, int nBefore, int newEdges, double timeSec) {
        IDFS *dfs = d.dfsStrategy;
        const std::vector<int> &rto = dfs->reverseTopologicalOrder;
        const std::vector<int> &topOrd = dfs->topOrd;
        int n = g.displayNumberVertices();

        // ---- current set, sorted keys
        std::vector<uint64_t> cur;
        cur.reserve(SB.size());
        for (const std::pair<int,int> &sb : SB) {
            if (sb.first == 0 && sb.second == 1) continue;
            cur.push_back(key(sb.first, sb.second));
        }
        std::sort(cur.begin(), cur.end());
        cur.erase(std::unique(cur.begin(), cur.end()), cur.end());
        // only superbubbles with an inside vertex (range entrance .. exit longer than the two)
        size_t all = cur.size();
        cur.erase(std::remove_if(cur.begin(), cur.end(),
                                 [&](uint64_t k) { return rto[second(k)] - rto[first(k)] <= 1; }), cur.end());
        long trivial = (long)(all - cur.size());

        // ---- nesting of the current set (sweep over the order)
        if ((int)exitOf.size() < n) { exitOf.resize(n, -1); parentOf.resize(n, -1); coverOf.resize(n, -1); }
        for (int s : entrances) { exitOf[s] = -1; parentOf[s] = -1; }
        entrances.clear();
        for (uint64_t k : cur) { int s = first(k); exitOf[s] = second(k); entrances.push_back(s); }
        std::vector<int> open;
        for (int p = 0; p < n; p++) {
            int x = topOrd[p];
            // close the superbubbles that ended before p; a superbubble starting
            // at p is not inside one ending at p (chain neighbours share p)
            while (!open.empty() && rto[exitOf[open.back()]] < p) open.pop_back();
            if (exitOf[x] != -1) {
                while (!open.empty() && rto[exitOf[open.back()]] <= p) open.pop_back();
                parentOf[x] = open.empty() ? -1 : open.back();
                open.push_back(x);
            }
            coverOf[x] = open.empty() ? -1 : open.back();   // innermost superbubble covering position p
        }
        // new vertices per position prefix (for "born" vs "promoted")
        newPrefix.assign(n + 1, 0);
        for (int p = 0; p < n; p++) newPrefix[p + 1] = newPrefix[p] + (topOrd[p] >= nBefore && nBefore > 0 ? 1 : 0);

        // ---- differences
        std::vector<uint64_t> born, died;
        std::set_difference(cur.begin(), cur.end(), prev.begin(), prev.end(), std::back_inserter(born));
        std::set_difference(prev.begin(), prev.end(), cur.begin(), cur.end(), std::back_inserter(died));

        long nBorn = 0, nPromoted = 0, nDiedMerged = 0;
        bornStamp.clear();
        for (uint64_t k : born) {
            int s = first(k), t = second(k);
            int inner = rto[t] - rto[s] - 1;
            bool hasNew = newPrefix[rto[t] + 1] - newPrefix[rto[s]] > 0;
            const char *ev = (step == 0) ? "base" : (hasNew ? "born" : "promoted");
            if (step > 0) { if (hasNew) nBorn++; else nPromoted++; }
            innerAtBirth[k] = inner;
            bornStamp[s] = 1;
            events << step << '\t' << genome << '\t' << ev << '\t' << g.toOriginal(s) << '\t' << g.toOriginal(t)
                   << '\t' << inner << '\t' << name(g, parentOf[s]) << '\n';
        }
        for (uint64_t k : died) {
            int s = first(k), t = second(k);
            int c = container(s, t, rto);
            if (c != -1 && bornStamp.count(c)) nDiedMerged++;
            std::unordered_map<uint64_t,int>::iterator it = innerAtBirth.find(k);
            int inner = (it == innerAtBirth.end()) ? -1 : it->second;
            if (it != innerAtBirth.end()) innerAtBirth.erase(it);
            events << step << '\t' << genome << '\t' << "died" << '\t' << g.toOriginal(s) << '\t' << g.toOriginal(t)
                   << '\t' << inner << '\t' << name(g, c) << '\n';
        }
        summary << step << '\t' << genome << '\t' << newEdges << '\t' << (nBefore > 0 ? n - nBefore : n - 2)
                << '\t' << (step == 0 ? (long)born.size() : nBorn) << '\t' << nPromoted << '\t' << died.size()
                << '\t' << nDiedMerged << '\t' << cur.size() << '\t' << timeSec << '\t' << trivial << '\n';
        prev.swap(cur);
    }

private:
    std::ofstream events, summary;
    std::vector<uint64_t> prev;                      // superbubbles after the previous step (sorted keys)
    std::unordered_map<uint64_t,int> innerAtBirth;   // inner vertex count when created
    std::vector<int> exitOf, parentOf, coverOf, entrances, newPrefix;
    std::unordered_map<int,int> bornStamp;           // entrances of the superbubbles born in this step

    static uint64_t key(int s, int t) { return ((uint64_t)(uint32_t)s << 32) | (uint32_t)t; }
    static int first(uint64_t k) { return (int)(k >> 32); }
    static int second(uint64_t k) { return (int)(k & 0xffffffffu); }

    // innermost current superbubble containing the positions of s and t
    int container(int s, int t, const std::vector<int> &rto) const {
        int a = std::min(rto[s], rto[t]), b = std::max(rto[s], rto[t]);
        int x = coverOf[s];
        while (x != -1 && !(rto[x] <= a && rto[exitOf[x]] >= b)) x = parentOf[x];
        return x;
    }

    std::string name(GraphContext &g, int entrance) const {
        if (entrance == -1) return "-";
        return std::to_string(g.toOriginal(entrance)) + "-" + std::to_string(g.toOriginal(exitOf[entrance]));
    }
};
#endif
