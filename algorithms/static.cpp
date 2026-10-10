#include <sys/resource.h>
#include "superbubbleAlgorithms.h" // Assuming this contains your class definitions
#include <fstream>
#include <chrono>
#include <string>
#include <sstream>

// LOG(...) is debug output: enabled only with `make ... DEBUG=1`
// (was: enabled unless BENCHMARK)
#ifdef DEBUG_OUTPUT
#define LOG(x) x
#else
#define LOG(x)
#endif

class StaticSuperbubble: public AlgorithmForSuperbubbles {
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

    void inputAndCompute(std::string output_path, std::string mode) {

        #ifndef BENCHMARK
        std::ofstream timeFile(output_path);
        #endif

        // --- Variables for parsing ---
        int num_base_edges, num_nodes, x, y;
        #ifndef BENCHMARK
        double minITime, maxOTime,insertTime, superbubbleTime, totalTime, currentinsertTime;
        #endif
        std::string line;
        std::string graph_name;
        int update_counter = 1;

        // --- Skip initial comment lines ---
        while (std::getline(std::cin, line) && (line.empty() || line[0] == '#')) {
            // This loop reads and discards all comment lines at the top of the file.
            //std::cout << line << std::endl;
        }
        //std::cout << line << std::endl;
        // --- 1. Process the Base Graph ---
        graph_name = line;
        std::cin >> num_base_edges >> num_nodes;

        // Initialize graph objects
        AdjacencyListGraph adl(num_base_edges); 
        GraphContext g(&adl);
        SimpleDFS simdfs(&g);
        DFSCONTEXT d(&simdfs);
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
        d.updateEdge(&g,0,1);
        for (int i = 0; i < num_base_edges; ++i) {
            std::cin >> x >> y;
            LOG(std::cout << x << " "<< y << std::endl;)
            d.updateEdge(&g, x, y); 
        }
        #ifdef DEBUG_OUTPUT  // debug output: compiled only with `make ... DEBUG=1` (was #ifndef BENCHMARK)
        std::cout << "All edges updated succesfully" << std::endl;
        #endif
        //g.displayOutgoingGraph();
        // Now, perform the first full calculation on the complete base graph
        #ifndef BENCHMARK
        auto start = std::chrono::high_resolution_clock::now();
        #endif
        #ifdef DEBUG_OUTPUT  // debug output: compiled only with `make ... DEBUG=1` (was #ifndef BENCHMARK)
        g.displayOutgoingGraph();
        #endif
        d.performDFS(&g);
        //std::cout << "DFS performed successfully" << std::endl;
        d.createMinI(&g);
        //std::cout << "minI done successfully" << std::endl;
        d.createMaxO(&g);
        //std::cout << "maxO done successfully" << std:: endl;
        #ifdef DEBUG_OUTPUT  // debug output: compiled only with `make ... DEBUG=1` (was #ifndef BENCHMARK)
        d.displayTopologicalOrder();
        d.displayReverseTopologicalOrder();
        d.displayMinI();
        d.displayMaxO();
        #endif
        superbubble.computeSB(&d, &g);
        //std::cout << "superbubble function done successfully" << std::endl;
        #ifdef DEBUG_OUTPUT  // debug output: compiled only with `make ... DEBUG=1` (was #ifndef BENCHMARK)
        superbubble.printSuperbubbles(&g);
        #endif
        #ifndef BENCHMARK
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        #endif
        
        
        


        #ifndef BENCHMARK
        timeFile << g.displayNumberEdges() << " " << g.displayNumberVertices() << "\n";
        timeFile << "BaseGraph_Build " << elapsed.count() << " " << superbubble.numberOfSuperbubbles() << "\n";

        std::cout << "base graph processed" << std::endl;
        #endif
        if(mode == "perEdge") {
            #ifndef BENCHMARK
            timeFile << "Edge_Number Haplotype_Name Insert_Time DFS_Time TopologicalOrderTime MinI_Time MaxO_Time SB_Compute_Time Total_Time SB_Count\n";
            #endif
        }
        else {
            #ifndef BENCHMARK
            timeFile << "Update_ID Haplotype_Name Insert_Time DFS_Time TopologicalOrderTime MinI_Time MaxO_Time SB_Compute_Time Total_Time SB_Count\n";
            #endif
        }
        // --- 2. Process Haplotype Additions Incrementally ---
        while (std::getline(std::cin, line)) {
            if (line.empty() || line.find("# HAPLOTYPE:") == std::string::npos) {
                continue;
            }
             
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
            #ifndef BENCHMARK
            minITime = 0.0;
            maxOTime = 0.0;
            insertTime = 0.0;
            currentinsertTime = 0.0;
            superbubbleTime = 0.0;
            totalTime = 0.0;
            #endif
            std::string clean_hap_name = hap_name_header.substr(hap_name_header.find(":") + 2);
            //start = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < hap_edge_count; ++i) {
                 std::cin >> x >> y;
                 #ifdef DEBUG_OUTPUT  // debug output: compiled only with `make ... DEBUG=1` (was #ifndef BENCHMARK)
                 std::cout << x << " " << y << std::endl;
                 #endif
                 #ifndef BENCHMARK
                 start = std::chrono::high_resolution_clock::now();
                 #endif
                 superbubble.insertE(&d,&g,x,y);
                 #ifdef DEBUG_OUTPUT  // debug output: compiled only with `make ... DEBUG=1` (was #ifndef BENCHMARK)
                 g.displayOutgoingGraph();
                 #endif
                 // use pointer to superbubble tree to flag affected superbubbles and all ancestors upto LCA
                 #ifndef BENCHMARK
                 end = std::chrono::high_resolution_clock::now();
                 elapsed = end - start;
                 currentinsertTime = elapsed.count();
                 insertTime += elapsed.count();
                 #endif
                 if (mode == "perEdge") {
                     
                     d.performDFS(&g);
                     #ifdef DEBUG_OUTPUT  // debug output: compiled only with `make ... DEBUG=1` (was #ifndef BENCHMARK)
                     d.displayTopologicalOrder();
                     d.displayReverseTopologicalOrder();
                     #endif
                     
                     

                
                     #ifndef BENCHMARK
                     start = std::chrono::high_resolution_clock::now();
                     #endif
                     d.createMinI(&g);
                     #ifndef BENCHMARK
                     end = std::chrono::high_resolution_clock::now();
                     elapsed = end - start;
                     minITime = elapsed.count();
                     #ifdef DEBUG_OUTPUT  // debug output: compiled only with `make ... DEBUG=1`
                     d.displayMinI();
                     #endif

                     start = std::chrono::high_resolution_clock::now();
                     #endif
                     d.createMaxO(&g);

                     #ifndef BENCHMARK
                     end = std::chrono::high_resolution_clock::now();
                     elapsed = end - start;
                     maxOTime = elapsed.count();
                     #ifdef DEBUG_OUTPUT  // debug output: compiled only with `make ... DEBUG=1`
                     d.displayMaxO();
                     #endif

                     start = std::chrono::high_resolution_clock::now();
                     #endif
                     superbubble.computeSB(&d,&g);
                     #ifndef BENCHMARK
                     end = std::chrono::high_resolution_clock::now();
                     elapsed = end - start;
                     superbubbleTime = elapsed.count();
                     #endif
                     
                     #ifdef DEBUG_OUTPUT  // debug output: compiled only with `make ... DEBUG=1` (was #ifndef BENCHMARK)
                     superbubble.printSuperbubbles(&g);
                     #endif
                     #ifndef BENCHMARK
                     totalTime = currentinsertTime + simdfs.latestTimes.dfsTime + simdfs.latestTimes.topologicalOrderTime + minITime + maxOTime + superbubbleTime;
                     timeFile << i + 1 << " " << clean_hap_name << " " << currentinsertTime << " " << simdfs.latestTimes.dfsTime << " " << simdfs.latestTimes.topologicalOrderTime << " " << minITime << " " << maxOTime << " " << superbubbleTime << " " << totalTime << " " << superbubble.numberOfSuperbubbles() << "\n";
                     currentinsertTime = 0.0;
                     simdfs.latestTimes.topologicalOrderTime = 0.0;
                     simdfs.latestTimes.dfsTime = 0.0;
                     #endif
                 }
                 

                 

            }
            
            //timeFile << g.displayNumberEdges() << " " << g.displayNumberVertices() << "\n";
            
            if (mode == "perBatch") {
                
                d.performDFS(&g);
                
                #ifndef BENCHMARK
                start = std::chrono::high_resolution_clock::now();
                #endif
                d.createMinI(&g); 
                #ifndef BENCHMARK
                end = std::chrono::high_resolution_clock::now();
                elapsed = end - start;
                minITime = elapsed.count();
                #endif
                
                
                #ifndef BENCHMARK
                start = std::chrono::high_resolution_clock::now();
                #endif
                d.createMaxO(&g);
                #ifndef BENCHMARK
                end = std::chrono::high_resolution_clock::now();
                elapsed = end - start;
                maxOTime = elapsed.count();
            

                start = std::chrono::high_resolution_clock::now();
                #endif
                superbubble.computeSB(&d,&g);
                #ifndef BENCHMARK
                end = std::chrono::high_resolution_clock::now();
                elapsed = end - start;
                superbubbleTime = elapsed.count();
                totalTime = insertTime + simdfs.latestTimes.dfsTime + simdfs.latestTimes.topologicalOrderTime + minITime + maxOTime + superbubbleTime;
                
                timeFile << update_counter << " " << clean_hap_name << " " << insertTime << " " << simdfs.latestTimes.dfsTime << " " << simdfs.latestTimes.topologicalOrderTime << " " << minITime << " " << maxOTime << " " << superbubbleTime << " " << totalTime << " " << superbubble.numberOfSuperbubbles() << "\n";
                simdfs.latestTimes.dfsTime = 0.0;
                simdfs.latestTimes.topologicalOrderTime = 0.0;
                #endif

            }
            #ifndef BENCHMARK
            timeFile << g.displayNumberEdges() << " " << g.displayNumberVertices() << "\n";
            #endif
            //timeFile << update_counter << " " << clean_hap_name << " " << elapsed.count() << " superbubble function" << "\n";
            //std::string clean_hap_name = hap_name_header.substr(hap_name_header.find(":") + 2);
            /*
            if (outFile.is_open()) {
                outFile << update_counter << " " << clean_hap_name << "\n";
                for (const pii& elem : superbubble.SB) {
                    outFile << elem.first << " " << elem.second << "\n";
                }
                std::cout << "List written to file successfully at: "<< std::endl;
                
            } else {
                std::cerr << "Error opening file at: " << std::endl; 
            }
            */
            /*
            totalTime = insertTime + dfsTime + minITime + maxOTime + superbubbleTime;
            timeFile << g.displayNumberEdges() << " " << g.displayNumberVertices() << "\n";
            timeFile << update_counter << " " << clean_hap_name << " " << insertTime << " " << dfsTime << " " << minITime << " " << maxOTime << " " << superbubbleTime << " " << totalTime << " " << superbubble.numberOfSuperbubbles() << "\n";
            */
            update_counter++;
        }
        //outFile.close();
        // --- Finalization ---
        #ifndef BENCHMARK
        timeFile << "PeakMemory(KB): " << getMemoryUsageInKB() << "\n";
        timeFile.close();
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
        std::cerr << "Usage: " << argv[0] << " <output_path> <perBatch|perEdge>" << std::endl;
        return 1; // Exit with an error code
    }
    StaticSuperbubble simple;
    simple.inputAndCompute(argv[1], argv[2]);
    return 0; 
}