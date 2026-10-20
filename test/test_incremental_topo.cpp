// Randomized correctness test + benchmark for incremental_topo.hpp
//   g++ -O2 -std=c++17 -DINCTOPO_DEBUG -I lib/DFS test/test_incremental_topo.cpp -o t && ./t
#include "incremental_topo.hpp"

#include <chrono>
#include <cstdio>
#include <numeric>
#include <random>

using namespace inctopo;

// Random batches drawn from a hidden permutation pi (edge a->b only if pi[a] < pi[b]),
// so every batch keeps the graph acyclic but generally violates the current order.
static bool correctness(uint64_t seed) {
    std::mt19937_64 rng(seed);
    long long batches = 0, edges = 0;
    for (int round = 0; round < 4000; ++round) {
        const int n = 2 + (int)(rng() % 60);
        std::vector<int> pi(n);
        std::iota(pi.begin(), pi.end(), 0);
        std::shuffle(pi.begin(), pi.end(), rng);

        IncrementalTopoOrder g;
        if (rng() % 3 == 0) {                         // sometimes start from a static build
            std::vector<Edge> init;
            for (int i = 0; i < n; ++i) {
                int a = (int)(rng() % n), b = (int)(rng() % n);
                if (a == b) continue;
                if (pi[a] > pi[b]) std::swap(a, b);
                init.push_back({a, b});
            }
            g.build(n, init);
            if (!g.verify()) { std::printf("build failed (round %d)\n", round); return false; }
        } else {
            g.reset(n);
        }

        const int nb = 1 + (int)(rng() % 20);
        for (int b = 0; b < nb; ++b) {
            if (rng() % 10 == 0) {                    // occasionally grow the graph
                Vertex v = g.addVertex();
                pi.push_back((int)(rng() % (pi.size() + 1)) - 1);  // any rank works; ties broken below
                (void)v;
                // re-rank pi to stay a strict order
                std::vector<int> idx(pi.size());
                std::iota(idx.begin(), idx.end(), 0);
                std::stable_sort(idx.begin(), idx.end(), [&](int x, int y) { return pi[x] < pi[y]; });
                for (size_t r = 0; r < idx.size(); ++r) pi[idx[r]] = (int)r;
            }
            const int nv = g.size();
            const int k = 1 + (int)(rng() % 12);
            std::vector<Edge> batch;
            for (int i = 0; i < k; ++i) {
                int a = (int)(rng() % nv), c = (int)(rng() % nv);
                if (a == c) continue;
                if (pi[a] > pi[c]) std::swap(a, c);
                batch.push_back({a, c});
            }
            g.insert(batch);
            ++batches; edges += (long long)batch.size();
            if (!g.verify()) {
                std::printf("FAIL: invalid order after batch (round %d, batch %d)\n", round, b);
                return false;
            }
        }
    }
    std::printf("correctness: %lld batches, %lld edges, all orders valid\n", batches, edges);
    return true;
}

static void benchmark(int n, long long m, int batchSize, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::vector<int> pi(n);
    std::iota(pi.begin(), pi.end(), 0);
    std::shuffle(pi.begin(), pi.end(), rng);

    std::vector<Edge> all;
    all.reserve(m);
    while ((long long)all.size() < m) {
        int a = (int)(rng() % n), b = (int)(rng() % n);
        if (a == b) continue;
        if (pi[a] > pi[b]) std::swap(a, b);
        all.push_back({a, b});
    }

    IncrementalTopoOrder g(n);
    long long window = 0, visited = 0;
    auto t0 = std::chrono::steady_clock::now();
    for (long long i = 0; i < m; i += batchSize) {
        const size_t k = (size_t)std::min<long long>(batchSize, m - i);
        BatchStats s = g.insert(all.data() + i, k);
        if (s.moved()) window += s.windowHi - s.windowLo + 1;
        visited += s.dfsVisited;
    }
    double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("bench n=%d m=%lld batch=%d: %.3f s (%.2f M edges/s), avg window %.0f, avg DFS visits %.0f, valid=%s\n",
                n, m, batchSize, sec, m / sec / 1e6,
                (double)window / ((m + batchSize - 1) / batchSize),
                (double)visited / ((m + batchSize - 1) / batchSize),
                g.verify() ? "yes" : "NO");
}

int main() {
    if (!correctness(12345)) return 1;
    if (!correctness(987654321)) return 1;
    // Sparse random DAGs, edges arriving in random order (adversarial-ish for the order).
    benchmark(100000, 300000, 1, 1);
    benchmark(100000, 300000, 100, 2);
    benchmark(100000, 300000, 10000, 3);
    benchmark(1000000, 3000000, 10000, 4);
    return 0;
}
