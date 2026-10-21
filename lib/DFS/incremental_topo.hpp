// incremental_topo.hpp
// Batch-incremental topological ordering of a DAG (DFS-based, bounded forward search).
// Header-only, C++17, no dependencies.
//
// Precondition: every inserted batch keeps the graph acyclic (no cycle checks
// in release builds; define INCTOPO_DEBUG to get asserts).
//
// State
//   ord_[v] : position of v in the topological order (0..n-1)
//   inv_[p] : vertex at position p (inverse of ord_, always kept in sync)
//   Invariant: for every edge (u,v), ord_[u] < ord_[v].
//
// insert(batch)
//   1. Append all edges to the adjacency lists (the searches must see them).
//      Only violating edges (ord[u] > ord[v]) need work; they are queued per
//      source vertex and define the window [lo, hi] = [min ord[v], max ord[u]].
//   2. Walk the window top-down (hi -> lo). For each vertex x not yet emitted,
//      run an iterative DFS from each queued target of x, restricted to
//      ord < ord[x]. Reached vertices are marked (shared across the batch) and
//      emitted in postorder; then x itself is emitted.
//   3. The emitted sequence is the new top-down order of the window:
//      out[k] gets position hi - k. Positions outside the window don't move.
//
// Cost per batch: O(k + |window| + edges scanned by the searches).
// Per-batch flags are epoch-stamped, so nothing is cleared between batches.
//
// Added for the superbubble pipeline (marked [added] / [changed] below):
//   assign()             start from a given order (the DFS order of the base
//                        graph) instead of build()'s Kahn order
//   addVertexAt()        place a new vertex where a DFS order needs it
//   assignWithNewVertices()  new order with many new vertices in one pass
//   moveBlock()          move a contiguous block of positions (a subtree)
//   insertAndRepair()    reorder for new edges and for edges already in the
//                        graph that point backwards, in one pass
//   movedInLastInsert()  which vertices the last reorder moved
// insert() is split into appending the edges and reorder(); same behaviour.
#pragma once

#include <algorithm>
#include <cstdint>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>
#ifdef INCTOPO_DEBUG
#include <cassert>
#define INCTOPO_ASSERT(x) assert(x)
#else
#define INCTOPO_ASSERT(x) ((void)0)
#endif

namespace inctopo {

using Vertex = int32_t;
using Edge = std::pair<Vertex, Vertex>;

struct BatchStats {
    int32_t windowLo = -1;    // first relabelled position (-1: nothing moved)
    int32_t windowHi = -1;    // last relabelled position (inclusive)
    int64_t dfsVisited = 0;   // vertices reached by the forward searches
    bool moved() const { return windowHi >= 0; }
};

class IncrementalTopoOrder {
public:
    IncrementalTopoOrder() = default;
    explicit IncrementalTopoOrder(Vertex n) { reset(n); }

    // n isolated vertices, identity order.
    void reset(Vertex n) {
        n = std::max<Vertex>(n, 0);
        adj_.assign(n, {});
        ord_.resize(n);
        inv_.resize(n);
        for (Vertex i = 0; i < n; ++i) ord_[i] = inv_[i] = i;
        mark_.assign(n, 0);
        pendStamp_.assign(n, 0);
        pendHead_.assign(n, -1);
        epoch_ = 0;
        numEdges_ = 0;
    }

    // Static build from an edge list (Kahn's algorithm, O(n + m)). Graph must be a DAG.
    void build(Vertex n, const std::vector<Edge>& edges) {
        reset(n);
        std::vector<int32_t> indeg(n, 0);
        for (const auto& e : edges) {
            adj_[e.first].push_back(e.second);
            ++indeg[e.second];
        }
        std::vector<Vertex> q;
        q.reserve(n);
        for (Vertex v = 0; v < n; ++v) if (indeg[v] == 0) q.push_back(v);
        for (size_t h = 0; h < q.size(); ++h)
            for (Vertex w : adj_[q[h]]) if (--indeg[w] == 0) q.push_back(w);
        INCTOPO_ASSERT((Vertex)q.size() == n && "build(): input is not a DAG");
        for (Vertex p = 0; p < (Vertex)q.size(); ++p) { inv_[p] = q[p]; ord_[q[p]] = p; }
        numEdges_ = edges.size();
    }

    // [added] Start from a given topological order: order[p] = vertex at
    // position p (e.g. a DFS order). Must be a valid order for the edges.
    void assign(Vertex n, const std::vector<Edge>& edges, const std::vector<Vertex>& order) {
        reset(n);
        for (const auto& e : edges) adj_[e.first].push_back(e.second);
        for (Vertex p = 0; p < n; ++p) { inv_[p] = order[p]; ord_[order[p]] = p; }
        numEdges_ = edges.size();
        INCTOPO_ASSERT(verify() && "assign(): order is not a topological order");
    }

    // New isolated vertex, appended at the end of the order (always valid).
    Vertex addVertex() {
        const Vertex v = size();
        adj_.emplace_back();
        ord_.push_back(v);
        inv_.push_back(v);
        mark_.push_back(0);
        pendStamp_.push_back(0);
        pendHead_.push_back(-1);
        return v;
    }

    // [added] New isolated vertex at position p; the vertices at positions
    // p.. move one position up (their relative order is unchanged). O(n - p).
    Vertex addVertexAt(int32_t p) {
        const Vertex v = size();
        INCTOPO_ASSERT(p >= 0 && p <= v);
        adj_.emplace_back();
        mark_.push_back(0);
        pendStamp_.push_back(0);
        pendHead_.push_back(-1);
        ord_.push_back(0);
        inv_.push_back(0);
        for (int32_t q = v; q > p; --q) { inv_[q] = inv_[q - 1]; ord_[inv_[q]] = q; }
        inv_[p] = v;
        ord_[v] = p;
        return v;
    }

    // [added] Adds the vertices size() .. order.size()-1 (no edges yet) and
    // takes 'order' as the new order: order[p] = vertex at position p; the old
    // vertices must keep their relative order. One pass, O(order.size()) -
    // for placing many new vertices at once instead of one addVertexAt each.
    void assignWithNewVertices(const std::vector<Vertex>& order) {
        const Vertex newN = (Vertex)order.size();
        INCTOPO_ASSERT(newN >= size());
        adj_.resize(newN);
        mark_.resize(newN, 0);
        pendStamp_.resize(newN, 0);
        pendHead_.resize(newN, -1);
        ord_.resize(newN);
        inv_.assign(order.begin(), order.end());
        for (Vertex p = 0; p < newN; ++p) ord_[inv_[p]] = p;
    }

    BatchStats insert(const std::vector<Edge>& batch) { return insert(batch.data(), batch.size()); }
    BatchStats insertEdge(Vertex u, Vertex v) { const Edge e{u, v}; return insert(&e, 1); }

    // [changed] split into appending the edges and reorder(); same behaviour
    BatchStats insert(const Edge* E, size_t k) {
        for (size_t i = 0; i < k; ++i) {
            INCTOPO_ASSERT(E[i].first >= 0 && E[i].first < size() && E[i].second >= 0 &&
                           E[i].second < size() && E[i].first != E[i].second);
            adj_[E[i].first].push_back(E[i].second);
        }
        numEdges_ += k;
        return reorder(E, k);
    }

    // [added] Appends newEdges like insert(), then repairs the order for the
    // new edges AND for 'existing' (edges already in the graph that point
    // backwards, e.g. after moveBlock) in one reorder pass.
    BatchStats insertAndRepair(const std::vector<Edge>& newEdges, const std::vector<Edge>& existing) {
        for (const Edge& e : newEdges) adj_[e.first].push_back(e.second);
        numEdges_ += newEdges.size();
        if (existing.empty()) return reorder(newEdges.data(), newEdges.size());
        all_.assign(newEdges.begin(), newEdges.end());
        all_.insert(all_.end(), existing.begin(), existing.end());
        return reorder(all_.data(), all_.size());
    }

    // [added] Moves the vertices at positions s..e (a contiguous block) to
    // positions to..to+(e-s), to <= s; the vertices at to..s-1 move behind the
    // block, keeping their order. O(e - to). Positions outside to..e unchanged.
    void moveBlock(int32_t s, int32_t e, int32_t to) {
        INCTOPO_ASSERT(0 <= to && to <= s && s <= e && e < size());
        std::rotate(inv_.begin() + to, inv_.begin() + s, inv_.begin() + e + 1);
        for (int32_t p = to; p <= e; ++p) ord_[inv_[p]] = p;
    }

private:
    // The reorder of insert(): the edges are already in the adjacency lists.
    BatchStats reorder(const Edge* E, size_t k) {
        BatchStats st;
        if (k == 0) return st;

        nextEpoch();
        const uint32_t ep = epoch_;
        pendNext_.resize(k);

        int32_t lo = std::numeric_limits<int32_t>::max(), hi = -1;
        for (size_t i = 0; i < k; ++i) {
            const Vertex u = E[i].first, v = E[i].second;
            if (ord_[u] > ord_[v]) {                          // violating edge: queue on u
                if (pendStamp_[u] != ep) { pendStamp_[u] = ep; pendHead_[u] = -1; }
                pendNext_[i] = pendHead_[u];
                pendHead_[u] = (int32_t)i;
                lo = std::min(lo, ord_[v]);
                hi = std::max(hi, ord_[u]);
            }
        }
        if (hi < 0) return st;                                 // every edge already consistent

        out_.clear();
        out_.reserve((size_t)(hi - lo + 1));
        for (int32_t p = hi; p >= lo; --p) {
            const Vertex x = inv_[p];
            if (mark_[x] == ep) continue;                      // already emitted by a search
            if (pendStamp_[x] == ep) {
                const int32_t bound = ord_[x];
                for (int32_t e = pendHead_[x]; e != -1; e = pendNext_[e])
                    st.dfsVisited += forward(E[e].second, bound, ep);
            }
            out_.push_back(x);
        }
        INCTOPO_ASSERT(out_.size() == (size_t)(hi - lo + 1) && "batch created a cycle");

        int32_t pos = hi;
        for (Vertex x : out_) { ord_[x] = pos; inv_[pos] = x; --pos; }
        st.windowLo = lo;
        st.windowHi = hi;
        return st;
    }

public:
    // --- queries -----------------------------------------------------------
    Vertex size() const { return (Vertex)ord_.size(); }
    size_t numEdges() const { return numEdges_; }
    int32_t position(Vertex v) const { return ord_[v]; }
    Vertex atPosition(int32_t p) const { return inv_[p]; }
    const std::vector<Vertex>& order() const { return inv_; }        // vertices in topological order
    const std::vector<int32_t>& positions() const { return ord_; }   // vertex -> position
    const std::vector<Vertex>& successors(Vertex v) const { return adj_[v]; }
    // Necessary condition for "u reaches v".
    bool mayReach(Vertex u, Vertex v) const { return ord_[u] < ord_[v]; }
    // [added] Was v reached by a forward search of the last insert(), i.e.
    // moved after the source of a violating edge? Only these vertices change
    // their position relative to the others; all other vertices of the
    // window keep their relative order.
    bool movedInLastInsert(Vertex v) const { return mark_[v] == epoch_; }

    void reserveVertices(size_t n) {
        adj_.reserve(n); ord_.reserve(n); inv_.reserve(n);
        mark_.reserve(n); pendStamp_.reserve(n); pendHead_.reserve(n);
    }

    // Full O(n + m) self-check (tests / debug).
    bool verify() const {
        const Vertex n = size();
        if ((Vertex)inv_.size() != n) return false;
        for (Vertex v = 0; v < n; ++v) {
            if (ord_[v] < 0 || ord_[v] >= n || inv_[ord_[v]] != v) return false;
            for (Vertex w : adj_[v]) if (ord_[v] >= ord_[w]) return false;
        }
        return true;
    }

private:
    struct Frame { Vertex v; uint32_t i; };

    void nextEpoch() {
        if (++epoch_ == 0) {   // wrap-around: clear stamps once every 2^32 batches
            std::fill(mark_.begin(), mark_.end(), 0);
            std::fill(pendStamp_.begin(), pendStamp_.end(), 0);
            epoch_ = 1;
        }
    }

    // Iterative DFS from v restricted to ord < bound; emits reached vertices in
    // postorder into out_. Returns the number of newly reached vertices.
    int64_t forward(Vertex v, int32_t bound, uint32_t ep) {
        if (mark_[v] == ep || ord_[v] >= bound) return 0;
        int64_t cnt = 1;
        mark_[v] = ep;
        stack_.clear();
        stack_.push_back({v, 0});
        while (!stack_.empty()) {
            Frame& f = stack_.back();
            const std::vector<Vertex>& a = adj_[f.v];
            if (f.i < a.size()) {
                const Vertex w = a[f.i++];
                if (mark_[w] == ep || ord_[w] >= bound) continue;
                mark_[w] = ep;
                ++cnt;
                stack_.push_back({w, 0});   // f may dangle after this; it is not used again
            } else {
                out_.push_back(f.v);
                stack_.pop_back();
            }
        }
        return cnt;
    }

    std::vector<std::vector<Vertex>> adj_;
    std::vector<int32_t> ord_, inv_;
    std::vector<uint32_t> mark_, pendStamp_;
    std::vector<int32_t> pendHead_, pendNext_;
    std::vector<Vertex> out_;
    std::vector<Edge> all_;                    // [added] insertAndRepair's combined edges
    std::vector<Frame> stack_;
    uint32_t epoch_ = 0;
    size_t numEdges_ = 0;
};

}  // namespace inctopo
