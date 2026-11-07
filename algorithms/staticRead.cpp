#include <sys/resource.h>
#include "superbubbleAlgorithms.h"
#include <fstream>
#include <chrono>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <memory>
#include <unordered_set>
#include <unistd.h>

/*
  staticRead: the static algorithm (static.cpp) run the way the external tools
  (BubbleGun, BubbleFinder) are run. After every update (an edge in perEdge,
  a haplotype in perBatch) it

    1. writes the current graph to a temporary GFA, exactly as
       scripts/superbubbleComparator.py writes it for the tools
       (Preprocessing/graphUpdater.py): the lines of the base GFA, an S line
       "S <name> *" for every segment the updates added, and an L line for
       every edge the updates added, in the order they were added
       (not timed, as for the tools);
    2. reads that GFA from disk into a new adjacency list (Read_Time: reading
       and parsing the file and inserting every edge, with the artificial
       source and sink);
    3. builds every data structure from scratch and finds the superbubbles
       as static does: DFS, topological order, minI, maxO, superbubbles.

  So nothing is kept from one update to the next except the update list
  itself, and the time is split like the tools' time: Read_Time (reading the
  graph) and Search_Time (everything after it, the work of static).

  Usage (run by scripts/superbubbleRunner.py -a staticRead):
      test/staticRead <time file> <perEdge|perBatch> <base GFA> [every] < numeric_edgelist.txt
  <base GFA> is reduced_graph.gfa next to the edge list. every (perEdge only,
  default 1): estimator as in static.cpp, the whole read + computation only
  after every <every>-th edge and after the last edge of every haplotype.

  Segment sX of the GFA is vertex X+1 of the edge list (0 and 1 are the
  artificial source and sink), as everywhere in the pipeline.
*/

// LOG(...) is debug output: enabled only with `make ... DEBUG=1`
#ifdef DEBUG_OUTPUT
#define LOG(x) x
#else
#define LOG(x)
#endif

// a graph and all data structures of static, built from scratch for every update
struct StaticBuild {
    AdjacencyListGraph adl;
    GraphContext g;
    SimpleDFS simdfs;
    DFSCONTEXT d;
    SuperbubbleCalculation superbubble;
    StaticBuild() : adl(0), g(&adl), simdfs(&g), d(&simdfs) {
        simdfs.clearOnInsert = false;       // one clear() before the DFS, as static perBatch
    }
};

// the times of one update (s)
struct StepTimes {
    double read = 0.0, clear = 0.0, dfs = 0.0, topo = 0.0, minI = 0.0, maxO = 0.0, sb = 0.0;
    double search() const { return dfs + topo + minI + maxO + sb; }
    double total() const { return read + clear + search(); }
};

class StaticReadSuperbubble: public AlgorithmForSuperbubbles {
public:
    // Peak resident memory (KB), as in static.cpp: VmHWM / VmRSS of this program
    long getCurrentMemoryInKB() {
        std::ifstream status("/proc/self/status");
        std::string line;
        while (std::getline(status, line)) {
            if (line.compare(0, 6, "VmRSS:") == 0) return std::stol(line.substr(6));
        }
        return -1;
    }
    long getMemoryUsageInKB() {
        std::ifstream status("/proc/self/status");
        std::string line;
        while (std::getline(status, line)) {
            if (line.compare(0, 6, "VmHWM:") == 0) return std::stol(line.substr(6));
        }
        return -1;
    }

    std::string baseGfa, tempGfa;
    std::vector<char> inBase;                    // vertex id -> segment has an S line in the base GFA
    std::vector<int> newNodes;                   // vertex ids added by the updates, in order
    std::unordered_set<int> newNodeSet;
    std::vector<std::pair<int,int>> addedEdges;  // edges added by the updates, in order
    std::unordered_set<long long> addedEdgeSet;

    // "sX" -> X+1; -1 if the name is not of that form
    static int vertexOf(const std::string &name) {
        if (name.size() < 2 || name[0] != 's') return -1;
        char *end = nullptr;
        long x = std::strtol(name.c_str() + 1, &end, 10);
        if (*end != '\0' || x < 0) return -1;
        return (int)(x + 1);
    }

    // the segments of the base GFA (once, not timed)
    void readBaseSegments() {
        std::ifstream in(baseGfa);
        if (!in.is_open()) {
            std::cerr << "staticRead: cannot open the base GFA " << baseGfa << std::endl;
            std::exit(1);
        }
        std::string line;
        while (std::getline(in, line)) {
            if (line.size() < 2 || line[0] != 'S' || line[1] != '\t') continue;
            size_t t = line.find('\t', 2);
            int v = vertexOf(line.substr(2, t == std::string::npos ? std::string::npos : t - 2));
            if (v < 0) {
                std::cerr << "staticRead: segment name not of the form s<number>: " << line.substr(0, 80) << std::endl;
                std::exit(1);
            }
            if ((int)inBase.size() <= v) inBase.resize(v + 1, 0);
            inBase[v] = 1;
        }
    }

    // an edge of an update (edge-list ids); as graphUpdater.py: duplicates are added once
    void addEdge(int x, int y) {
        long long key = ((long long)x << 32) | (unsigned int)y;
        if (!addedEdgeSet.insert(key).second) return;
        addedEdges.push_back(std::make_pair(x, y));
        for (int v : {x, y}) {
            bool base = v < (int)inBase.size() && inBase[v];
            if (!base && newNodeSet.insert(v).second) newNodes.push_back(v);
        }
    }

    // the temporary GFA of the current graph (not timed), as graphUpdater.py writes it
    void writeTempGfa() {
        std::ofstream out(tempGfa, std::ios::binary | std::ios::trunc);
        {
            std::ifstream in(baseGfa, std::ios::binary);
            out << in.rdbuf();
        }
        for (int v : newNodes) out << "S\ts" << v - 1 << "\t*\n";
        for (const auto &e : addedEdges) out << "L\ts" << e.first - 1 << "\t+\ts" << e.second - 1 << "\t+\t0M\n";
        out.close();
        if (!out) {
            std::cerr << "staticRead: cannot write " << tempGfa << std::endl;
            std::exit(1);
        }
    }

    // reads the temporary GFA into a new graph and finds its superbubbles as static does
    std::unique_ptr<StaticBuild> readAndCompute(StepTimes &t) {
        auto start = std::chrono::high_resolution_clock::now();
        std::unique_ptr<StaticBuild> b(new StaticBuild());
        b->d.updateEdge(&b->g, 0, 1);
        std::ifstream in(tempGfa);
        std::string line;
        while (std::getline(in, line)) {
            if (line.size() < 2 || line[0] != 'L' || line[1] != '\t') continue;     // S lines: read, not needed
            // L <from> <+> <to> <+> <overlap>
            size_t a = 2, a2 = line.find('\t', a);
            size_t o1 = a2 + 1, b1 = line.find('\t', o1) + 1, b2 = line.find('\t', b1);
            size_t o2 = b2 + 1;
            if (a2 == std::string::npos || b2 == std::string::npos || line[o1] != '+' || line[o2] != '+') {
                std::cerr << "staticRead: L line not of the form L s<a> + s<b> + : " << line.substr(0, 80) << std::endl;
                std::exit(1);
            }
            int x = vertexOf(line.substr(a, a2 - a)), y = vertexOf(line.substr(b1, b2 - b1));
            if (x < 0 || y < 0) {
                std::cerr << "staticRead: segment name not of the form s<number>: " << line.substr(0, 80) << std::endl;
                std::exit(1);
            }
            b->d.updateEdge(&b->g, x, y);
        }
        auto end = std::chrono::high_resolution_clock::now();
        t.read = std::chrono::duration<double>(end - start).count();

        // the data structures of static, from scratch
        start = std::chrono::high_resolution_clock::now();
        b->simdfs.clear();
        end = std::chrono::high_resolution_clock::now();
        t.clear = std::chrono::duration<double>(end - start).count();
        b->d.performDFS(&b->g);
        t.dfs = b->simdfs.latestTimes.dfsTime;
        t.topo = b->simdfs.latestTimes.topologicalOrderTime;
        start = std::chrono::high_resolution_clock::now();
        b->d.createMinI(&b->g);
        end = std::chrono::high_resolution_clock::now();
        t.minI = std::chrono::duration<double>(end - start).count();
        start = std::chrono::high_resolution_clock::now();
        b->d.createMaxO(&b->g);
        end = std::chrono::high_resolution_clock::now();
        t.maxO = std::chrono::duration<double>(end - start).count();
        start = std::chrono::high_resolution_clock::now();
        b->superbubble.computeSB(&b->d, &b->g);
        end = std::chrono::high_resolution_clock::now();
        t.sb = std::chrono::duration<double>(end - start).count();
        #ifdef DEBUG_OUTPUT
        b->superbubble.printSuperbubbles(&b->g);
        #endif
        return b;
    }

    // one row of the time file: the columns of static, then Read_Time and Search_Time
    // (UpdateDS_Time = reading the GFA into the adjacency list + clearing the arrays)
    static void writeRow(std::ofstream &f, long long id, const std::string &hap, const StepTimes &t, int sbCount) {
        f << id << " " << hap << " " << t.read + t.clear << " " << t.dfs << " " << t.topo << " " << t.minI << " "
          << t.maxO << " " << t.sb << " " << t.total() << " " << sbCount << " " << t.read << " " << t.search();
    }

    void inputAndCompute(std::string output_path, std::string mode, std::string base_gfa, int every = 1) {
        long baselineMemoryKB = getCurrentMemoryInKB();   // before the input is read
        baseGfa = base_gfa;
        char tmpl[] = "/tmp/staticRead_XXXXXX";
        int fd = mkstemp(tmpl);
        if (fd < 0) {
            std::cerr << "staticRead: cannot create a temporary file" << std::endl;
            std::exit(1);
        }
        close(fd);
        tempGfa = tmpl;
        readBaseSegments();

        #ifndef BENCHMARK
        std::ofstream timeFile(output_path);
        if (!timeFile.is_open()) {
            std::cerr << "Error opening time file!" << std::endl;
            return;
        }
        #endif

        int num_base_edges, num_nodes, x, y;
        std::string line;
        int update_counter = 1;
        bool sampling = (mode == "perEdge" && every > 1);
        long long globalStep = 0, evenSamples = 0, otherSamples = 0;
        double evenTime = 0.0;

        // --- Skip initial comment lines; graph name ---
        while (std::getline(std::cin, line) && (line.empty() || line[0] == '#')) {}
        std::cin >> num_base_edges >> num_nodes;
        // the base edges are in the base GFA; read past them
        for (int i = 0; i < num_base_edges; ++i) std::cin >> x >> y;

        // --- 1. The base graph ---
        StepTimes t;
        writeTempGfa();
        std::unique_ptr<StaticBuild> b = readAndCompute(t);
        #ifndef BENCHMARK
        timeFile << b->g.displayNumberEdges() << " " << b->g.displayNumberVertices() << "\n";
        timeFile << "BaseGraph_Build " << t.total() << " " << b->superbubble.numberOfSuperbubbles() << "\n";
        timeFile << (mode == "perEdge" ? "Edge_Number" : "Update_ID")
                 << " Haplotype_Name UpdateDS_Time DFS_Time TopologicalOrderTime MinI_Time MaxO_Time SB_Compute_Time"
                    " Total_Time SB_Count Read_Time Search_Time" << (sampling ? " Global_Step Sample" : "") << "\n";
        std::cout << "base graph processed" << std::endl;
        #endif

        // --- 2. The updates ---
        while (std::getline(std::cin, line)) {
            if (line.empty() || line.find("# HAPLOTYPE:") == std::string::npos) continue;
            size_t colon_pos = line.find(":");
            size_t name_start_pos = line.find_first_not_of(" \t", colon_pos + 1);
            if (name_start_pos == std::string::npos || line[name_start_pos] == '#') continue;
            std::string clean_hap_name = line.substr(line.find(":") + 2);
            int hap_edge_count = 0;
            if (std::cin.peek() >= '0' && std::cin.peek() <= '9') std::cin >> hap_edge_count;
            if (hap_edge_count == 0) {
                #ifndef BENCHMARK
                timeFile << update_counter << " " << clean_hap_name << " " << 0.0 << " "
                         << b->superbubble.numberOfSuperbubbles() << "\n";
                #endif
                update_counter++;
                std::string dummy;
                std::getline(std::cin, dummy);
                continue;
            }
            for (int i = 0; i < hap_edge_count; ++i) {
                std::cin >> x >> y;
                addEdge(x, y);
                if (mode != "perEdge") continue;
                globalStep++;
                bool evenSample = (globalStep % every == 0);
                if (sampling && !evenSample && i != hap_edge_count - 1) continue;
                writeTempGfa();
                b.reset();                               // the previous graph is freed, not timed
                b = readAndCompute(t);
                #ifndef BENCHMARK
                writeRow(timeFile, i + 1, clean_hap_name, t, b->superbubble.numberOfSuperbubbles());
                if (sampling) {
                    timeFile << " " << globalStep << " " << (evenSample ? "E" : "L");
                    if (evenSample) { evenSamples++; evenTime += t.total(); }
                    else otherSamples++;
                }
                timeFile << "\n";
                #endif
            }
            if (mode == "perBatch") {
                writeTempGfa();
                b.reset();
                b = readAndCompute(t);
                #ifndef BENCHMARK
                writeRow(timeFile, update_counter, clean_hap_name, t, b->superbubble.numberOfSuperbubbles());
                timeFile << "\n";
                #endif
            }
            #ifndef BENCHMARK
            timeFile << b->g.displayNumberEdges() << " " << b->g.displayNumberVertices() << "\n";
            #endif
            update_counter++;
        }

        // --- Finalization ---
        std::remove(tempGfa.c_str());
        #ifndef BENCHMARK
        if (sampling) {
            double mean = evenSamples > 0 ? evenTime / evenSamples : 0.0;
            timeFile << "Updates: " << globalStep << "\n";
            timeFile << "Every: " << every << "\n";
            timeFile << "SampledUpdates: " << evenSamples << "\n";
            timeFile << "RunsAfterLastEdgeOfGenome: " << otherSamples << "\n";
            timeFile << "MeanUpdateTime: " << mean << "\n";
            timeFile << "TotalTimeEstimated: " << mean * globalStep << "\n";
        }
        timeFile << "BaselineMemory(KB): " << baselineMemoryKB << "\n";
        timeFile << "PeakMemory(KB): " << getMemoryUsageInKB() << "\n";
        timeFile.close();
        // superbubbles of the final graph -> staticRead_superbubbles.txt
        writeFinalSuperbubbles(output_path, b->superbubble, b->g);
        #endif
    }
};

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <output_path> <perBatch|perEdge> <base GFA> [every] < numeric_edgelist.txt"
                  << std::endl;
        std::cerr << "  every (perEdge only): estimator, read and compute only after every <every>-th edge" << std::endl;
        return 1;
    }
    int every = (argc >= 5) ? std::atoi(argv[4]) : 1;
    if (every < 1) every = 1;
    StaticReadSuperbubble simple;
    simple.inputAndCompute(argv[1], argv[2], argv[3], every);
    return 0;
}
