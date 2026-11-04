#include <sys/resource.h>
#include"superbubbleAlgorithms.h"
#include<fstream>
#include<string>
#include<chrono>
// LOG(...) is debug output: enabled only with `make ... DEBUG=1`
// (was: enabled unless BENCHMARK)
#ifdef DEBUG_OUTPUT
#define LOG(x) x
#else
#define LOG(x)
#endif
/*
  incBSTBatchH2: the batch algorithm (incBSTBatch: batch DFS, minI, maxO per
  genome) with the superbubble heuristic per batch restricted to the ranges of
  the groups (SuperbubbleTree::updateBatchRanges), so P = root does not force a
  full computation as in incBSTBatchH1. perBatch only.
*/
class IncBSTBatchH2Superbubble: public AlgorithmForSuperbubbles {
public:  
    #ifdef VERIFY_SB
    // make ... VERIFY=1: after every genome, the superbubbles must equal
    // those of the full routine (SuperbubbleCalculation::computeSB)
    long verifyChecked = 0, verifyFailed = 0;
    void verifySuperbubbles(SuperbubbleCalculationBatchH1 &sb, DFSCONTEXT &d, GraphContext &g, const std::string &genome) {
        SuperbubbleCalculation ref;
        ref.computeSB(&d, &g);
        std::vector<pii> a = sb.superbubbles();
        std::vector<pii> b(ref.SB.begin(), ref.SB.end());
        std::sort(a.begin(), a.end());
        std::sort(b.begin(), b.end());
        verifyChecked++;
        if (a == b) return;
        verifyFailed++;
        if (verifyFailed <= 5) {
            std::cout << "VERIFY: genome " << verifyChecked << " (" << genome << ") differs from the full routine" << std::endl;
            int shown = 0;
            for (auto &e : a) if (!std::binary_search(b.begin(), b.end(), e) && shown++ < 10)
                std::cout << "  only incremental: <" << g.toOriginal(e.first) << " " << g.toOriginal(e.second) << ">" << std::endl;
            shown = 0;
            for (auto &e : b) if (!std::binary_search(a.begin(), a.end(), e) && shown++ < 10)
                std::cout << "  only full:        <" << g.toOriginal(e.first) << " " << g.toOriginal(e.second) << ">" << std::endl;
        }
        // continue from the correct state
        sb.computeSB(&d, &g);
    }
    #endif

    // Peak resident memory (KB) of this program: VmHWM from /proc/self/status.
    // Not getrusage(RUSAGE_SELF).ru_maxrss: on Linux a program started by
    // another process inherits that process's peak, so when started from the
    // Python runners it reported the runner's memory instead of its own.
    // resident memory now (VmRSS): read at the start, before the input, it is the
    // program's own memory without the graph (baseline of PeakMemory)
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
            if (line.compare(0, 6, "VmHWM:") == 0) {
                return std::stol(line.substr(6));   // "VmHWM:     4960 kB"
            }
        }
        return -1;
    }   
    void inputAndCompute(std::string output_path, std::string mode) {
        #ifndef BENCHMARK
        long baselineMemoryKB = getCurrentMemoryInKB();   // before the input is read
        #endif
        // --- File path configuration ---
        //std::ofstream outFile("/home/dell/Research/Superbubble/Experiment/IncSuperbubble/output/Human/HaplotypeTestingData/outFile/SuperbubbleIncrementalAVLDFSOutput.txt");
        #ifndef BENCHMARK
        std::ofstream timeFile(output_path);
        #endif
        // --- Variables for parsing ---
        int num_base_edges, num_nodes, x, y;
        #ifndef BENCHMARK
        double superbubbleTime, currentsuperbubbleTime, totalTime;
        #endif
        std::string line;
        std::string graph_name;
        int update_counter = 1;

        // --- Skip initial comment lines ---
        while (std::getline(std::cin, line) && (line.empty() || line[0] == '#')) {
            // This loop reads and discards all comment lines at the top of the file.
            #ifdef DEBUG_OUTPUT  // debug output: compiled only with `make ... DEBUG=1` (was #ifndef BENCHMARK)
            std::cout << line << std::endl;
            #endif
        }
        //std::cout << line << std::endl;
        // --- 1. Process the Base Graph ---
        graph_name = line;
        std::cin >> num_base_edges >> num_nodes;
        LOG(std::cout << num_base_edges << " " << num_nodes << std::endl;)
        // Initialize graph objects
        // <-- CORRECTED LINE: Initialize with the number of edges, as per your class design.
        AdjacencyListGraph adl(num_base_edges); 
        GraphContext g(&adl);
        IncrementalDFSTREEBatch incdfs(&g);
        DFSCONTEXT d(&incdfs);
        SuperbubbleCalculationBatchH1 superbubble(true);     // by ranges
        LOG(std::cout << "creating of objects classes" << std::endl;)
        // Write initial header to time file
        
        #ifndef BENCHMARK
        if (!timeFile.is_open()) {
            std::cerr << "Error opening time file!" << std::endl;
            return;
        }
        std::cout << num_base_edges << std::endl;
        #endif

        // BATCH INSERT all edges from the base graph first.
        for (int i = 0; i < num_base_edges; ++i) {
            std::cin >> x >> y;
            LOG(std::cout << x << " "<< y << std::endl;)
            d.updateEdge(&g, x, y); 
        }

        // Now, perform the first full calculation on the complete base graph
        #ifndef BENCHMARK
        auto start = std::chrono::high_resolution_clock::now();
        #endif
        d.performDFS(&g);   // static DFS of the base graph
        d.createMinI(&g);
        d.createMaxO(&g);
        superbubble.computeSB(&d, &g);
        #ifndef BENCHMARK
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        #endif

        //std::cout << "DFS performed successfully" << std::endl;
        //std::cout << "minI done successfully" << std::endl;
        //std::cout << "maxO done successfully" << std:: endl;
        //d.displayTopologicalOrder();
        //d.displayReverseTopologicalOrder();
        //d.displayMinI();
        //d.displayMaxO();
        #ifdef DEBUG_OUTPUT  // debug output: compiled only with `make ... DEBUG=1` (was #ifndef BENCHMARK)
        superbubble.printSuperbubbles(&g);
        #endif
        

        // Record timing for the entire base graph construction
        #ifndef BENCHMARK
        timeFile << g.displayNumberEdges() << " " << g.displayNumberVertices() << "\n";
        timeFile << "BaseGraph_Build " << elapsed.count() << " " << superbubble.numberOfSuperbubbles() << "\n";
        
        std::cout << "base graph processed" << std::endl;
        #endif
        #ifndef BENCHMARK
        // per genome: what the batch did (NewVertices SourcesPlaced BackwardEdges
        // Rounds Adoptions, see IncrementalDFSTREEBatch::BatchCounts) instead of
        // the per-edge cases of the other algorithms
        timeFile << "Update_ID Haplotype_Name UpdateDS_Time DFS_Time TopologicalOrderTime MinI_Time MaxO_Time SB_Compute_Time Total_Time SB_Count NewVertices SourcesPlaced BackwardEdges Rounds Adoptions AffectedVertices ExploredVertices TotalVertices TotalEdges TotalAffectedEdges TotalExploredEdges\n";
        #endif
        //timeFile << "Update_ID Haplotype_Name UpdateDS_Time DFS_Time MinI_Time MaxO_Time SB_Compute_Time Total_Time SB_Count Case1 Case2 Case3 Case4 Case5\n";
        // --- 2. Process Haplotype Additions Incrementally ---
        while (std::getline(std::cin, line)) {
            if (line.empty() || line.find("# HAPLOTYPE:") == std::string::npos) {
                continue;
            }
             
            // --- NEW, SMARTER CHECK ---
            // Find the position of the colon
            size_t colon_pos = line.find(":");
            // Find the first non-space character after the colon
            size_t name_start_pos = line.find_first_not_of(" \t", colon_pos + 1);

            // If there's nothing after the colon, or if the first character is a '#',
            // then it's a separator line, so we skip it.
            if (name_start_pos == std::string::npos || line[name_start_pos] == '#') {
              continue;
            }
            // --- END OF NEW CHECK ---

            std::string hap_name_header = line;
            int hap_edge_count = 0;
            
            // Carefully read the next line to get the edge count.
            if (std::cin.peek() >= '0' && std::cin.peek() <= '9') {
                 std::cin >> hap_edge_count;
            }
            
            if (hap_edge_count == 0) {
                 std::string clean_hap_name = hap_name_header.substr(hap_name_header.find(":") + 2);
                 #ifndef BENCHMARK
                 timeFile << update_counter << " " << clean_hap_name << " " << 0.0 << " " << superbubble.numberOfSuperbubbles() << "\n";
                 #endif
                 update_counter++;
                 std::string dummy;
                 std::getline(std::cin, dummy);
                 continue;
            }
            std::string clean_hap_name = hap_name_header.substr(hap_name_header.find(":") + 2);
            #ifndef BENCHMARK
            std::cout << clean_hap_name << std::endl;

            superbubbleTime = 0.0;
            currentsuperbubbleTime = 0.0;
            totalTime = 0.0;
            incdfs.latestTimes.insertTime = 0.0;
            incdfs.latestTimes.dfsTime = 0.0;
            incdfs.latestTimes.topologicalOrderTime = 0.0;
            incdfs.latestTimes.minITime = 0.0;
            incdfs.latestTimes.maxOTime = 0.0;
            incdfs.vertexCount.countAffectedVertices = 0;
            incdfs.vertexCount.countExploredVertices = 0;
            incdfs.vertexCount.totalEdges = 0;
            incdfs.vertexCount.totalExploredEdges = 0;
            incdfs.vertexCount.totalAffectedEdges = 0;
            #endif
            // Start the timer for this entire haplotype update
            
            
             // BATCH INSERT all edges for this haplotype
            std::vector<std::pair<int,int>> edges;           // the genome's edges, applied as one batch
            edges.reserve(hap_edge_count);
            for (int i = 0; i < hap_edge_count; ++i) {
                std::cin >> x >> y;
                edges.push_back(std::make_pair(x, y));
            }
            
            if (mode == "perBatch") {
                // the whole haplotype as one batch
                superbubble.insertBatch(&d,&g,edges);
                #ifndef BENCHMARK
                start = std::chrono::high_resolution_clock::now();
                #endif
                superbubble.updateSB(&d, &g);
                #ifndef BENCHMARK
                end = std::chrono::high_resolution_clock::now();
                // + reading the groups before the batch (in insertBatch)
                superbubbleTime = std::chrono::duration<double>(end - start).count() + superbubble.prepareTime;
                totalTime = incdfs.latestTimes.insertTime + incdfs.latestTimes.dfsTime + incdfs.latestTimes.topologicalOrderTime + incdfs.latestTimes.minITime + incdfs.latestTimes.maxOTime + superbubbleTime;
                timeFile << update_counter << " " << clean_hap_name << " " << incdfs.latestTimes.insertTime << " " << incdfs.latestTimes.dfsTime << " " << incdfs.latestTimes.topologicalOrderTime << " " << incdfs.latestTimes.minITime << " " << incdfs.latestTimes.maxOTime << " " << superbubbleTime << " " << totalTime << " " << superbubble.numberOfSuperbubbles() << " " << incdfs.counts.currentNewVertices << " " << incdfs.counts.currentSourcesPlaced << " " << incdfs.counts.currentBackwardEdges << " " << incdfs.counts.currentRounds << " " << incdfs.counts.currentAdoptions << " " << incdfs.vertexCount.countAffectedVertices << " " << incdfs.vertexCount.countExploredVertices << " " << incdfs.n << " " << incdfs.vertexCount.totalEdges << " " << incdfs.vertexCount.totalAffectedEdges << " " << incdfs.vertexCount.totalExploredEdges << "\n";
                #endif 
                #ifdef VERIFY_SB
                verifySuperbubbles(superbubble, d, g, clean_hap_name);
                #endif
            }
            //std::string clean_hap_name = hap_name_header.substr(hap_name_header.find(":") + 2);
            //std::cout << clean_hap_name << std::endl;
            /*
            if (outFile.is_open()) {
                outFile << update_counter << " " << clean_hap_name << "\n";
                for (const pii& elem : superbubble.SB) {
                    outFile << elem.first << " " << elem.second << "\n";
                }
                std::cout << "List written to file successfully at: "<< std::endl;
                
            } else {
                std::cerr << "Error opening file at: " << std::endl; 
            }*/
            //totalTime = incdfs.latestTimes.insertTime + incdfs.latestTimes.dfsTime + incdfs.latestTimes.minITime + incdfs.latestTimes.maxOTime + superbubbleTime;
            #ifndef BENCHMARK
            timeFile << g.displayNumberEdges() << " " << g.displayNumberVertices() << "\n";
            #endif
            //timeFile << update_counter << " " << clean_hap_name << " " << incdfs.latestTimes.insertTime << " " << incdfs.latestTimes.dfsTime << " " << incdfs.latestTimes.minITime << " " << incdfs.latestTimes.maxOTime << " " << superbubbleTime << " " << totalTime << " " << superbubble.numberOfSuperbubbles() << " " << incdfs.edgeTypes.case1 << " " << incdfs.edgeTypes.case2 << " " << incdfs.edgeTypes.case3 << " " << incdfs.edgeTypes.case4 << " " << incdfs.edgeTypes.case5 << "\n";
            update_counter++;
        }
        //outFile.close();
        // --- Finalization ---
        #ifndef BENCHMARK
        timeFile << "BaselineMemory(KB): " << baselineMemoryKB << "\n";
        timeFile << "PeakMemory(KB): " << getMemoryUsageInKB() << "\n";
        timeFile.close();
        // superbubbles of the final graph -> <ID>_superbubbles.txt
        writeFinalSuperbubbles(output_path, superbubble, g);
        std::cout << "batch algorithm: " << incdfs.counts.newVertices << " new vertices, " << incdfs.counts.sourcesPlaced << " sources placed under their first parent, "
                  << incdfs.counts.backwardEdges << " backward edges in " << incdfs.counts.rounds << " reorder rounds, " << incdfs.counts.adoptions << " adoptions ("
                  << incdfs.repairedEdges << " backward edges from them), " << incdfs.duplicateEdges << " duplicate edges" << std::endl;
        {
            const SuperbubbleTree::Stats &st = superbubble.tree.stats;
            std::cout << "superbubble routine: " << st.batches << " batches, " << st.batchesNoScan << " without a scan, "
                      << st.rescans << " scans (" << st.rescanVertices << " vertices scanned), "
                      << st.fullRecomputes << " full recomputations (" << st.fullSinkSource << " tail sink / head source, "
                      << st.fullRootP << " P = root, " << st.fullDangling << " dangling, " << st.fullNotClosed << " region not closed), "
                      << superbubble.duplicateEdges << " duplicate edges; "
                      << st.added << " superbubbles added, " << st.killed << " killed" << std::endl;
        }
        #endif
        #ifdef VERIFY_SB
        std::cout << "VERIFY: " << verifyChecked << " genomes checked, " << verifyFailed << " differ from a full recomputation" << std::endl;
        #endif
    }

};
int main(int argc, char* argv[]){
    // argc is the argument count, argv are the arguments
    // argv[0] is the program name
    // argv[1] will be our output path

    if (argc < 3) {
        // Print error to standard error
        std::cerr << "Error: No output file path provided." << std::endl;
        std::cerr << "Usage: " << argv[0] << " <output_file_path> perBatch" << std::endl;
        return 1; // Exit with an error code
    }
    // the batch algorithm applies one genome at a time: perBatch only
    if (std::string(argv[2]) != "perBatch") {
        std::cerr << "Error: This implementation only supports batch updates." << std::endl;
        std::cerr << "Please run with mode = perBatch." << std::endl;
        return 1;
    }
	IncBSTBatchH2Superbubble incremental;
    incremental.inputAndCompute(argv[1], argv[2]);
    return 0;
}
