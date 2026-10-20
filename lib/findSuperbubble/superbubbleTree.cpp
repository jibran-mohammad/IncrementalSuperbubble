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
    }
    return id;
}

void SuperbubbleTree::kill(int id) {
    assert(id != ROOT && nalive[id]);
    nalive[id] = 0;
    if (byEntrance[ns[id]] == id) byEntrance[ns[id]] = -1;
    if (byExit[nt[id]] == id) byExit[nt[id]] = -1;
    --aliveCount; ++stats.killed;
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
// ---------------------------------------------------------------------------
void SuperbubbleTree::rescan(const SBView& V, int owner) {
    int L = V.pos[ns[owner]], R = V.pos[nt[owner]];
    assert(L < R);

    // 1. kill the old superbubbles inside (every one has its entrance strictly inside)
    for (int p = L + 1; p < R; ++p) {
        int x = V.at[p];
        if (byEntrance[x] >= 0) kill(byEntrance[x]);
    }

    // 2. scan
    std::vector<std::pair<int,int>> pairs = scanRange(V, L, R);
    std::vector<int> exitAt(R - L + 1, -1);      // indexed by entrance position - L
    bool ownerFound = false;
    for (auto& pr : pairs) {
        if (pr.first == L && pr.second == R) { ownerFound = true; continue; }
        assert(pr.first > L);                    // otherwise owner was not a superbubble
        exitAt[pr.first - L] = pr.second;
    }
    assert(ownerFound);
    (void)ownerFound;

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
}

void SuperbubbleTree::build(const SBView& V) {
    ns.clear(); nt.clear(); npar.clear(); nmark.clear(); nalive.clear();
    aliveCount = 0;
    byEntrance.assign(V.n, -1);
    byExit.assign(V.n, -1);
    encl.assign(V.n, ROOT);
    newNode(0, 1, -1);                          // ROOT
    assert(V.pos[0] == 0 && V.pos[1] == V.n - 1);
    rescan(V, ROOT);
}

void SuperbubbleTree::rescanOrFull(const SBView& V, int P) {
    if (P == ROOT) { ++stats.fullRecomputes; ++stats.fullRootP; build(V); }
    else rescan(V, P);
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
