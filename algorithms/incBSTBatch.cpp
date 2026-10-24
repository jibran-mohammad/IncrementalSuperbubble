#include <sys/resource.h>
#include"superbubbleAlgorithms.h"
#include"../lib/findSuperbubble/superbubbleHistory.h"
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
class IncBSTBatchSuperbubble: public AlgorithmForSuperbubbles {
public:  
    // Peak resident memory (KB) of this program: VmHWM from /proc/self/status.
    // Not getrusage(RUSAGE_SELF).ru_maxrss: on Linux a program started by
    // another process inherits that process's peak, so when started from the
    // Python runners it reported the runner's memory instead of its own.
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
    // history: also write the superbubble history per genome (perBatch only,
    // normal build): <ID>_history_events.tsv, <ID>_history_summary.tsv next to
    // the time file (see lib/findSuperbubble/superbubbleHistory.h)
    void inputAndCompute(std::string output_path, std::string mode, bool history) {
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
        SuperbubbleCalculation superbubble;
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
        SuperbubbleHistory hist;
        if (history) {
            std::string prefix = output_path;
            const std::string suffix = "_time.txt";
            if (prefix.size() >= suffix.size() && prefix.compare(prefix.size() - suffix.size(), suffix.size(), suffix) == 0)
                prefix = prefix.substr(0, prefix.size() - suffix.size());
            if (!hist.open(prefix)) { std::cerr << "history: cannot write " << prefix << "_history_*.tsv" << std::endl; history = false; }
            else hist.record(0, "base", d, g, superbubble.SB, 0, num_base_edges, elapsed.count());
        }
        #endif
        #ifndef BENCHMARK
        timeFile << "Update_ID Haplotype_Name Insert_Time DFS_Time TopologicalOrderTime MinI_Time MaxO_Time SB_Compute_Time Total_Time SB_Count Case1 Case2 Case3 Case4 Case5 AffectedVertices ExploredVertices TotalVertices TotalEdges TotalAffectedEdges TotalExploredEdges\n";
        #endif
        //timeFile << "Update_ID Haplotype_Name Insert_Time DFS_Time MinI_Time MaxO_Time SB_Compute_Time Total_Time SB_Count Case1 Case2 Case3 Case4 Case5\n";
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
                 if (history) hist.record(update_counter, clean_hap_name, d, g, superbubble.SB, g.displayNumberVertices(), 0, 0.0);
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
            incdfs.edgeTypes.case1 = 0;
            incdfs.edgeTypes.case2 = 0;
            incdfs.edgeTypes.case3 = 0;
            incdfs.edgeTypes.case4 = 0;
            incdfs.edgeTypes.case5 = 0;
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
                #ifndef BENCHMARK
                int nBefore = g.displayNumberVertices();     // for the history: new vertices have ids >= nBefore
                #endif
                superbubble.insertBatch(&d,&g,edges);
                #ifndef BENCHMARK
                start = std::chrono::high_resolution_clock::now();
                #endif
                superbubble.computeSB(&d, &g);
                #ifndef BENCHMARK
                end = std::chrono::high_resolution_clock::now();
                superbubbleTime = std::chrono::duration<double>(end - start).count();
                totalTime = incdfs.latestTimes.insertTime + incdfs.latestTimes.dfsTime + incdfs.latestTimes.topologicalOrderTime + incdfs.latestTimes.minITime + incdfs.latestTimes.maxOTime + superbubbleTime;
                timeFile << update_counter << " " << clean_hap_name << " " << incdfs.latestTimes.insertTime << " " << incdfs.latestTimes.dfsTime << " " << incdfs.latestTimes.topologicalOrderTime << " " << incdfs.latestTimes.minITime << " " << incdfs.latestTimes.maxOTime << " " << superbubbleTime << " " << totalTime << " " << superbubble.numberOfSuperbubbles() << " " << incdfs.edgeTypes.case1 << " " << incdfs.edgeTypes.case2 << " " << incdfs.edgeTypes.case3 << " " << incdfs.edgeTypes.case4 << " " << incdfs.edgeTypes.case5 << " "<< incdfs.vertexCount.countAffectedVertices << " " << incdfs.vertexCount.countExploredVertices << " " << incdfs.n << " " << incdfs.vertexCount.totalEdges << " " << incdfs.vertexCount.totalAffectedEdges << " " << incdfs.vertexCount.totalExploredEdges << "\n";
                if (history) hist.record(update_counter, clean_hap_name, d, g, superbubble.SB, nBefore, hap_edge_count, totalTime);
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
        timeFile << "PeakMemory(KB): " << getMemoryUsageInKB() << "\n";
        timeFile.close();
        // superbubbles of the final graph -> <ID>_superbubbles.txt
        writeFinalSuperbubbles(output_path, superbubble, g);
        std::cout << "batch algorithm: " << incdfs.rounds << " reorder rounds, " << incdfs.sourcesPlaced << " sources placed under their first parent, " << incdfs.adoptions << " adoptions (case 5, v's only parent 0; "
                  << incdfs.repairedEdges << " backward edges repaired), " << incdfs.duplicateEdges << " duplicate edges" << std::endl;
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
        std::cerr << "Usage: " << argv[0] << " <output_file_path> perBatch [history]" << std::endl;
        return 1; // Exit with an error code
    }
    // the batch algorithm applies one genome at a time: perBatch only
    if (std::string(argv[2]) != "perBatch") {
        std::cerr << "Error: This implementation only supports batch updates." << std::endl;
        std::cerr << "Please run with mode = perBatch." << std::endl;
        return 1;
    }
	IncBSTBatchSuperbubble incremental;
    bool history = (argc >= 4 && std::string(argv[3]) == "history");
    incremental.inputAndCompute(argv[1], argv[2], history);
    return 0;
}
