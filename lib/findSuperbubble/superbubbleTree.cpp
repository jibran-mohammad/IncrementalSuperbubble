#include "superbubbleTree.h"
#include <algorithm>
#include <cassert>

namespace {
// pinned sentinels: vertex 1 behaves as if its earliest parent were 0,
// vertex 0 as if its latest child were 1. In the incremental DFS minI[0] and
// maxO[1] are INT_MAX (no real parent / child), so both are pinned for both
// sentinels (0 is only ever at position 0, 1 only at position n-1).
inline int minPos(const SBView& V, int x) { return (x == 0 || x == 1) ? 0 : V.pos[V.minI[x]]; }
inline int maxPos(const SBView& V, int x) { return (x == 0 || x == 1) ? V.n - 1 : V.pos[V.maxO[x]]; }
}

// ---------------------------------------------------------------------------
// The stack algorithm of SuperbubbleCalculation::computeSB, in positions and
// restricted to [L,R]. Each stack entry is (entrance position, max child position
// of the entries merged into it).
// ---------------------------------------------------------------------------
std::vector<std::pair<int,int>> SuperbubbleTree::scanRange(const SBView& V, int L, int R) const {
    struct E { int s, mx; };
    std::vector<E> S;
    std::vector<std::pair<int,int>> out;
    for (int p = L; p <= R; ++p) {
        int x = V.at[p];
        int mi = (p == L) ? L : minPos(V, x);       // L's parents lie outside the range
        // entries that start after x's earliest parent cannot be entrances any more
        while (!S.empty() && S.back().s > mi) {
            E top = S.back(); S.pop_back();
            if (!S.empty()) S.back().mx = std::max(S.back().mx, top.mx);
        }
        // the top entry closes exactly at x: smallest exit for that entrance
        if (!S.empty() && S.back().mx == p) {
            out.emplace_back(S.back().s, p);
            E top = S.back(); S.pop_back();
            if (!S.empty()) S.back().mx = std::max(S.back().mx, top.mx);
        }
        S.push_back({p, maxPos(V, x)});
    }
    return out;
}

// ---------------------------------------------------------------------------
// node bookkeeping
// ---------------------------------------------------------------------------
void SuperbubbleTree::ensureSize(int n) {
    if ((int)byEntrance.size() < n) {
        byEntrance.resize(n, -1);
        byExit.resize(n, -1);
        encl.resize(n, ROOT);       // a new vertex is a source or a sink: only the root contains it
    }
}

int SuperbubbleTree::newNode(int s, int t, int parent) {
    int id = (int)ns.size();
    ns.push_back(s); nt.push_back(t);
    npar.push_back(parent < 0 ? id : parent);
    nalive.push_back(1); nmark.push_back(0);
    if (id != ROOT) {
        assert(byEntrance[s] < 0 && byExit[t] < 0);   // one superbubble per entrance / exit
        byEntrance[s] = id; byExit[t] = id;
        ++aliveCount; ++stats.added;
        if (journal) journalAdded.emplace_back(s, t);
    }
    return id;
}

void SuperbubbleTree::kill(int id) {
    assert(id != ROOT && nalive[id]);
    nalive[id] = 0;
    if (byEntrance[ns[id]] == id) byEntrance[ns[id]] = -1;
    if (byExit[nt[id]] == id) byExit[nt[id]] = -1;
    --aliveCount; ++stats.killed;
    if (journal) journalKilled.emplace_back(ns[id], nt[id]);
}

// nearest alive node at or above id; dead nodes keep their parent pointers
int SuperbubbleTree::findAlive(int id) {
    int r = id;
    while (!nalive[r]) r = npar[r];
    while (id != r) { int nx = npar[id]; npar[id] = r; id = nx; }   // path compression
    return r;
}

int SuperbubbleTree::innerOut(int x) {
    if (byEntrance[x] >= 0) return byEntrance[x];
    return encl[x] = findAlive(encl[x]);
}

int SuperbubbleTree::innerIn(int x) {
    if (byExit[x] >= 0) return byExit[x];
    return encl[x] = findAlive(encl[x]);
}

// climbs from a and b alternately; cost ~ 2 * (dist(a,P) + dist(b,P))
int SuperbubbleTree::lca(int a, int b) {
    ++stamp;
    int x = a, y = b;
    while (true) {
        if (x >= 0) {
            if (nmark[x] == stamp) return x;
            nmark[x] = stamp;
            x = (x == ROOT) ? -1 : findAlive(npar[x]);
        }
        if (y >= 0) {
            if (nmark[y] == stamp) return y;
            nmark[y] = stamp;
            y = (y == ROOT) ? -1 : findAlive(npar[y]);
        }
    }
}

void SuperbubbleTree::killPath(int from, int to) {
    int x = from;
    while (x != to) {
        assert(x != ROOT);
        int up = findAlive(npar[x]);
        kill(x);
        x = up;
    }
}

// ---------------------------------------------------------------------------
// Replace everything strictly inside 'owner' by a fresh scan of its range.
// Requires: owner's interval is closed in the current graph (range lemma).
// Returns false, and changes nothing, if owner is no longer a superbubble
// (only possible in updateBatch; the per-edge cases guarantee it).
// ---------------------------------------------------------------------------
bool SuperbubbleTree::rescan(const SBView& V, int owner) {
    int L = V.pos[ns[owner]], R = V.pos[nt[owner]];
    if (L >= R) return false;

    // 1. scan
    std::vector<std::pair<int,int>> pairs = scanRange(V, L, R);
    std::vector<int> exitAt(R - L + 1, -1);      // indexed by entrance position - L
    bool ownerFound = false;
    for (auto& pr : pairs) {
        if (pr.first == L && pr.second == R) { ownerFound = true; continue; }
        if (pr.first <= L) return false;         // owner is not a superbubble
        exitAt[pr.first - L] = pr.second;
    }
    if (!ownerFound) return false;

    // 2. kill the old superbubbles inside (every one has its entrance strictly inside)
    for (int p = L + 1; p < R; ++p) {
        int x = V.at[p];
        if (byEntrance[x] >= 0) kill(byEntrance[x]);
    }

    // 3. sweep left to right: create nodes, set parents and encl
    std::vector<int> st;                         // open nodes, innermost on top
    for (int p = L; p <= R; ++p) {
        while (!st.empty() && V.pos[nt[st.back()]] <= p) st.pop_back();
        int x = V.at[p];
        int top = st.empty() ? owner : st.back();
        if (p > L && p < R) encl[x] = top;
        if (exitAt[p - L] >= 0) {
            int id = newNode(x, V.at[exitAt[p - L]], top);
            st.push_back(id);
        }
    }
    ++stats.rescans;
    stats.rescanVertices += R - L + 1;
    return true;
}

void SuperbubbleTree::build(const SBView& V) {
    if (journal)                                // the pool is dropped: every alive superbubble is removed
        for (int id = 1; id < (int)ns.size(); ++id)
            if (nalive[id]) journalKilled.emplace_back(ns[id], nt[id]);
    ns.clear(); nt.clear(); npar.clear(); nmark.clear(); nalive.clear();
    aliveCount = 0;
    byEntrance.assign(V.n, -1);
    byExit.assign(V.n, -1);
    encl.assign(V.n, ROOT);
    newNode(0, 1, -1);                          // ROOT
    assert(V.pos[0] == 0 && V.pos[1] == V.n - 1);
    bool ok = rescan(V, ROOT);
    assert(ok);
    (void)ok;
}

void SuperbubbleTree::rescanOrFull(const SBView& V, int P) {
    if (P == ROOT) { ++stats.fullRecomputes; ++stats.fullRootP; build(V); }
    else {
        bool ok = rescan(V, P);
        assert(ok);
        (void)ok;
    }
}

// ---------------------------------------------------------------------------
// The five cases. Proofs are in the conversation notes; summary:
//   1  u new          : v had a parent -> In(v) chain dies, nothing new
//                       v was a source -> only (u,v) is added
//   2  v new          : u had a child  -> Out(u) chain dies, nothing new
//                       u was a sink   -> only (u,v) is added
//   3  both new       : (u,v) is added
//   4  p(u) > p(v)    : u sink or v source -> full; otherwise delete both
//                       paths below P = LCA(innerOut(u), innerIn(v)), no scan
//   5  p(u) < p(v)    : u sink or v source -> full; otherwise kill v-side path,
//                       and if the u-side path is non-empty, rescan P
// ---------------------------------------------------------------------------
void SuperbubbleTree::update(const SBView& V, int caseNo, int u, int v,
                             bool uWasSink, bool vWasSource) {
    ensureSize(V.n);
    switch (caseNo) {
    case 1:
        if (vWasSource) newNode(u, v, ROOT);
        else killPath(innerIn(v), ROOT);
        ++stats.pathOnly;
        break;
    case 2:
        if (uWasSink) newNode(u, v, ROOT);
        else killPath(innerOut(u), ROOT);
        ++stats.pathOnly;
        break;
    case 3:
        newNode(u, v, ROOT);
        ++stats.pathOnly;
        break;
    case 4:
        if (uWasSink || vWasSource) { ++stats.fullRecomputes; ++stats.fullSinkSource; build(V); break; }
        {   // nothing new can appear: delete both paths below P, no scan
            int a = innerOut(u), b = innerIn(v), P = lca(a, b);
            killPath(a, P); killPath(b, P); ++stats.pathOnly;
        }
        break;
    case 5: {
        if (uWasSink || vWasSource) { ++stats.fullRecomputes; ++stats.fullSinkSource; build(V); break; }
        int a = innerOut(u), b = innerIn(v);
        int P = lca(a, b);
        killPath(b, P);                         // v-side: destroyed, never replaced
        if (a == P) { ++stats.pathOnly; break; } // u-side empty: nothing can be created
        rescanOrFull(V, P);                     // u-side destroyed and replaced inside P
        break;
    }
    default:
        assert(false && "unknown case");
    }
}

// ---------------------------------------------------------------------------
// A batch (one genome) at once, incBSTBatchH1. The groups are read from the
// graph before the batch (SuperbubbleCalculationBatchH1::insertBatch); the
// tree is still the one before the batch, V is the order after it.
//   closed group (old tails and old heads): a tail was a sink or a head a
//       source -> full (as 4b / 5c / 5d); P = LCA of innerOut(tails) and
//       innerIn(heads) on the old tree; P = root -> full (as 5a);
//       otherwise P is rescanned
//   dangling group (no old head or no old tail): edge by edge it is cases
//       1-3, so it must be a simple path of new vertices (otherwise full):
//       a tail that was not a sink loses every superbubble around it up to
//       the root, a head that was not a source the same (killPath, cases 1
//       and 2); consecutive path vertices are new superbubbles (case 3), and
//       so are (tail, front) / (back, head) if tailPair / headPair
//   then: a P killed by a dangling group -> full; only the outermost P's are
//   rescanned (a P inside another one is rebuilt by its rescan); a P that is
//   no longer a superbubble after the batch -> full
// ---------------------------------------------------------------------------
void SuperbubbleTree::updateBatch(const SBView& V, const std::vector<BatchGroup>& groups) {
    ensureSize(V.n);
    ++stats.batches;

    // 1. closed groups: their region P on the old tree
    std::vector<int> regions;
    bool dangling = false;
    for (const BatchGroup& gr : groups) {
        if (gr.tails.empty() || gr.heads.empty()) {
            if (!gr.isPath) { fullBuild(V, stats.fullDangling); return; }
            dangling = true;
            continue;
        }
        if (gr.tailWasSink || gr.headWasSource) { fullBuild(V, stats.fullSinkSource); return; }
        int P = innerOut(gr.tails[0]);
        for (size_t i = 1; i < gr.tails.size(); ++i) P = lca(P, innerOut(gr.tails[i]));
        for (int h : gr.heads) P = lca(P, innerIn(h));
        if (P == ROOT) { fullBuild(V, stats.fullRootP); return; }
        regions.push_back(P);
    }

    // 2. dangling groups: the superbubbles around their old vertex die (cases 1, 2)
    if (dangling) {
        for (const BatchGroup& gr : groups) {
            if (!gr.tails.empty() && !gr.heads.empty()) continue;
            for (int t : gr.tails) if (!gr.tailWasSink) killPath(innerOut(t), ROOT);
            for (int h : gr.heads) if (!gr.headWasSource) killPath(innerIn(h), ROOT);
        }
        for (int P : regions) if (!nalive[P]) { fullBuild(V, stats.fullDangling); return; }
    }

    // 3. rescan the outermost regions, from left to right in the new order
    std::sort(regions.begin(), regions.end(), [&](int a, int b) { return V.pos[ns[a]] < V.pos[ns[b]]; });
    regions.erase(std::unique(regions.begin(), regions.end()), regions.end());
    int end = -1;                               // exit position of the last rescanned region
    for (int P : regions) {
        if (V.pos[ns[P]] < end) continue;       // inside that region: already rebuilt
        if (!rescan(V, P)) { fullBuild(V, stats.fullNotClosed); return; }
        end = V.pos[nt[P]];
    }
    if (regions.empty()) ++stats.batchesNoScan;

    // 4. the superbubbles along the dangling paths (case 3; tailPair / headPair: case 2 / case 1)
    if (dangling) {
        for (const BatchGroup& gr : groups) {
            if (!gr.tails.empty() && !gr.heads.empty()) continue;
            const std::vector<int>& p = gr.path;
            if (gr.tailPair) addLeaf(gr.tails[0], p.front());
            for (size_t i = 0; i + 1 < p.size(); ++i) addLeaf(p[i], p[i + 1]);
            if (gr.headPair) addLeaf(p.back(), gr.heads[0]);
        }
    }
}

// ---------------------------------------------------------------------------
// A batch by ranges, incBSTBatchH2. As updateBatch, but a closed group does
// not rescan its whole region P: only the part of P between the group's ends.
//   a tail / head directly inside P counts with its own position, otherwise
//   with both ends of the child of P that holds it: the children outside stay
//   superbubbles, nothing in them was touched
//   ranges of the same P that overlap or touch are merged (a new superbubble
//   can span two touched children that share a vertex), so P = root needs no
//   full computation any more
//   full computation: a tail was a sink or a head a source (a superbubble can
//   appear anywhere), a dangling group that is not a simple path, a closed
//   group with a new vertex without a child or a parent (it opens the
//   superbubbles around the group, and a range scan does not see that), a
//   range that cuts a superbubble that stays (rescanRange returns false)
// The dangling groups kill first, so the P's are alive.
// ---------------------------------------------------------------------------
int SuperbubbleTree::childUnder(int x, int P) {
    int up = findAlive(npar[x]);
    while (up != P) { x = up; up = findAlive(npar[x]); }
    return x;
}

// Replace the superbubbles with their entrance in [L,R) by a fresh scan of that
// range, under 'owner' (alive, contains the range). Changes nothing and returns
// false if the range does not fit the tree.
bool SuperbubbleTree::rescanRange(const SBView& V, int owner, int L, int R) {
    if (L >= R) return true;
    std::vector<std::pair<int,int>> pairs = scanRange(V, L, R);

    // 1. the old superbubbles in the range must end in it
    for (int p = L; p < R; ++p) {
        int id = byEntrance[V.at[p]];
        if (id >= 0 && id != owner && V.pos[nt[id]] > R) return false;
    }
    for (int p = L + 1; p <= R; ++p) {
        int id = byExit[V.at[p]];
        if (id >= 0 && id != owner && V.pos[ns[id]] < L) return false;
    }
    // 2. the new ones must fit next to what stays
    std::vector<int> exitAt(R - L + 1, -1);      // indexed by entrance position - L
    for (auto& pr : pairs) {
        int s = V.at[pr.first], t = V.at[pr.second];
        if (s == ns[owner]) {
            if (t == nt[owner]) continue;        // the owner itself
            return false;                        // the owner is no longer a superbubble
        }
        if (t == nt[owner]) return false;
        int other = byExit[t];                   // stays if its entrance is outside [L,R)
        if (other >= 0 && (V.pos[ns[other]] < L || V.pos[ns[other]] >= R)) return false;
        exitAt[pr.first - L] = pr.second;
    }

    // 3. kill the old ones
    for (int p = L; p < R; ++p) {
        int id = byEntrance[V.at[p]];
        if (id >= 0 && id != owner) kill(id);
    }

    // 4. sweep left to right: create nodes, set parents and encl
    std::vector<int> st;                         // open nodes, innermost on top
    for (int p = L; p <= R; ++p) {
        while (!st.empty() && V.pos[nt[st.back()]] <= p) st.pop_back();
        int x = V.at[p];
        int top = st.empty() ? owner : st.back();
        if (p > L && p < R) encl[x] = top;
        if (exitAt[p - L] >= 0) {
            int id = newNode(x, V.at[exitAt[p - L]], top);
            st.push_back(id);
        }
    }
    ++stats.rescans;
    stats.rescanVertices += R - L + 1;
    return true;
}

void SuperbubbleTree::updateBatchRanges(const SBView& V, const std::vector<BatchGroup>& groups) {
    ensureSize(V.n);
    ++stats.batches;

    // 1. what needs a full computation
    bool dangling = false;
    for (const BatchGroup& gr : groups) {
        if (gr.tails.empty() || gr.heads.empty()) {
            if (!gr.isPath) { fullBuild(V, stats.fullDangling); return; }
            dangling = true;
        }
        else if (gr.tailWasSink || gr.headWasSource) { fullBuild(V, stats.fullSinkSource); return; }
        else if (gr.openEnd) { fullBuild(V, stats.fullDangling); return; }
    }

    // 2. dangling groups: the superbubbles around their old vertex die (cases 1, 2)
    if (dangling) {
        for (const BatchGroup& gr : groups) {
            if (!gr.tails.empty() && !gr.heads.empty()) continue;
            for (int t : gr.tails) if (!gr.tailWasSink) killPath(innerOut(t), ROOT);
            for (int h : gr.heads) if (!gr.headWasSource) killPath(innerIn(h), ROOT);
        }
    }

    // 3. closed groups: P on the tree and the range inside P, in the new order
    struct Range { int L, R, P; };
    std::vector<Range> ranges;
    for (const BatchGroup& gr : groups) {
        if (gr.tails.empty() || gr.heads.empty()) continue;
        int P = innerOut(gr.tails[0]);
        for (size_t i = 1; i < gr.tails.size(); ++i) P = lca(P, innerOut(gr.tails[i]));
        for (int h : gr.heads) P = lca(P, innerIn(h));
        Range r = { V.n, -1, P };
        auto cover = [&](int x, int inner) {
            r.L = std::min(r.L, V.pos[x]);
            r.R = std::max(r.R, V.pos[x]);
            if (inner == P) return;
            int c = childUnder(inner, P);
            r.L = std::min(r.L, V.pos[ns[c]]);
            r.R = std::max(r.R, V.pos[nt[c]]);
        };
        for (int t : gr.tails) cover(t, innerOut(t));
        for (int h : gr.heads) cover(h, innerIn(h));
        ranges.push_back(r);
    }

    // 4. merge from left to right and rescan
    std::sort(ranges.begin(), ranges.end(), [](const Range& a, const Range& b) {
        return a.L != b.L ? a.L < b.L : a.R > b.R;
    });
    std::vector<Range> merged;
    for (const Range& r : ranges) {
        if (!merged.empty()) {
            Range& m = merged.back();
            if (r.L < m.R || (r.L == m.R && r.P == m.P)) {
                if (r.P != m.P) m.P = lca(m.P, r.P);
                m.R = std::max(m.R, r.R);
                continue;
            }
        }
        merged.push_back(r);
    }
    for (const Range& m : merged) {
        if (!rescanRange(V, m.P, m.L, m.R)) { fullBuild(V, stats.fullNotClosed); return; }
    }
    if (merged.empty()) ++stats.batchesNoScan;

    // 5. the superbubbles along the dangling paths (case 3; tailPair / headPair: case 2 / case 1)
    if (dangling) {
        for (const BatchGroup& gr : groups) {
            if (!gr.tails.empty() && !gr.heads.empty()) continue;
            const std::vector<int>& p = gr.path;
            if (gr.tailPair) addLeaf(gr.tails[0], p.front());
            for (size_t i = 0; i + 1 < p.size(); ++i) addLeaf(p[i], p[i + 1]);
            if (gr.headPair) addLeaf(p.back(), gr.heads[0]);
        }
    }
}

void SuperbubbleTree::fullBuild(const SBView& V, long long& reason) {
    ++stats.fullRecomputes;
    ++reason;
    build(V);
}

// a superbubble without inside vertices, directly under the root
void SuperbubbleTree::addLeaf(int s, int t) {
    if (byEntrance[s] < 0 && byExit[t] < 0) newNode(s, t, ROOT);
}

std::vector<std::pair<int,int>> SuperbubbleTree::superbubbles() const {
    std::vector<std::pair<int,int>> out;
    for (int id = 1; id < (int)ns.size(); ++id)
        if (nalive[id]) out.emplace_back(ns[id], nt[id]);
    return out;
}

bool SuperbubbleTree::checkAgainstFullScan(const SBView& V) const {
    std::vector<std::pair<int,int>> fresh;
    for (auto& pr : scanRange(V, 0, V.n - 1))
        if (!(pr.first == 0 && pr.second == V.n - 1))
            fresh.emplace_back(V.at[pr.first], V.at[pr.second]);
    std::vector<std::pair<int,int>> mine = superbubbles();
    std::sort(fresh.begin(), fresh.end());
    std::sort(mine.begin(), mine.end());
    return fresh == mine;
}
