#include <sys/resource.h>
#include "superbubbleAlgorithms.h" // Assuming this contains your class definitions
#include <fstream>
#include <chrono>
#include <string>
#include <sstream>

class SimpleSuperbubbleH1: public AlgorithmForSuperbubbles {
public:
    long getMemoryUsageInKB() {
        struct rusage usage;
        getrusage(RUSAGE_SELF, &usage);
        return usage.ru_maxrss; // memory in kilobytes
    }

    void inputAndCompute(std::string output_path, std::string mode) {

        std::ofstream timeFile(output_path);

        // --- Variables for parsing ---
        int num_base_edges, num_nodes, x, y;
        double dfsTime, minITime, maxOTime,insertTime, superbubbleTime, totalTime, currentinsertTime;
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
        SuperbubbleCalculationH1 superbubble;
        std::cout << "creating of objects classes" << std::endl;
        // Write initial header to time file
        
        if (!timeFile.is_open()) {
            std::cerr << "Error opening time file!" << std::endl;
            return;
        }
        std::cout << num_base_edges << std::endl;
        
        // BATCH INSERT all edges from the base graph first.
        d.updateEdge(&g,0,1);
        for (int i = 0; i < num_base_edges; ++i) {
            std::cin >> x >> y;
            std::cout << x << " "<< y << std::endl;
            d.updateEdge(&g, x, y); 
        }

        // Now, perform the first full calculation on the complete base graph
        auto start = std::chrono::high_resolution_clock::now();
        d.performDFS(&g);
        d.createMinI(&g);
        d.createMaxO(&g);
        superbubble.computeSB(&d, &g);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        /*
        outFile << graph_name << "\n";
        if (outFile.is_open()) {
            for (const pii& elem : superbubble.SB) {
                outFile << elem.first << " " << elem.second << "\n";
            }
            std::cout << "List written to file successfully at: "<< std::endl;
        } else {
            std::cerr << "Error opening file at: " << std::endl; 
        }
        */
        // Record timing for the entire base graph construction

        timeFile << g.displayNumberEdges() << " " << g.displayNumberVertices() << "\n";
        timeFile << "BaseGraph_Build " << elapsed.count() << " " << superbubble.numberOfSuperbubbles() << "\n";

        std::cout << "base graph processed" << std::endl;
        if(mode == "perEdge") {
            timeFile << "Edge_Number Haplotype_Name Insert_Time DFS_Time MinI_Time MaxO_Time SB_Compute_Time Total_Time SB_Count Case1 Case2 Case3 Case4 Case5\n";
        }
        else {
            timeFile << "Update_ID Haplotype_Name Insert_Time DFS_Time MinI_Time MaxO_Time SB_Compute_Time Total_Time SB_Count Case1 Case2 Case3 Case4 Case5\n";
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
                 timeFile << update_counter << " " << clean_hap_name << " " << 0.0 << " " << superbubble.numberOfSuperbubbles() << "\n";
                 update_counter++;
                 std::string dummy;
                 std::getline(std::cin, dummy);
                 continue;
            }
            dfsTime = 0.0;
            minITime = 0.0;
            maxOTime = 0.0;
            insertTime = 0.0;
            currentinsertTime = 0.0;
            superbubbleTime = 0.0;
            totalTime = 0.0;
            std::string clean_hap_name = hap_name_header.substr(hap_name_header.find(":") + 2);
            //start = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < hap_edge_count; ++i) {
                 std::cin >> x >> y;
                 start = std::chrono::high_resolution_clock::now();
                 superbubble.insertE(&d,&g,x,y);
                 // use pointer to superbubble tree to flag affected superbubbles and all ancestors upto LCA
                 end = std::chrono::high_resolution_clock::now();
                 elapsed = end - start;
                 currentinsertTime = elapsed.count();
                 insertTime += elapsed.count();
                 
                 if (mode == "perEdge") {
                     start = std::chrono::high_resolution_clock::now();
                     d.performDFS(&g);
                     end = std::chrono::high_resolution_clock::now();
                     elapsed = end - start;
                     dfsTime = elapsed.count();

                     start = std::chrono::high_resolution_clock::now();
                     d.createMinI(&g);
                     end = std::chrono::high_resolution_clock::now();
                     elapsed = end - start;
                     minITime = elapsed.count();

                     start = std::chrono::high_resolution_clock::now();
                     d.createMaxO(&g);

                     end = std::chrono::high_resolution_clock::now();
                     elapsed = end - start;
                     maxOTime = elapsed.count();

                     start = std::chrono::high_resolution_clock::now();
                     superbubble.recomputeInside(&d,&g,x,y);
                     end = std::chrono::high_resolution_clock::now();
                     elapsed = end - start;
                     superbubbleTime = elapsed.count();
                     
                     totalTime = currentinsertTime + dfsTime + minITime + maxOTime + superbubbleTime;
                     timeFile << i + 1 << " " << clean_hap_name << " " << currentinsertTime << " " << dfsTime << " " << minITime << " " << maxOTime << " " << superbubbleTime << " " << totalTime << " " << superbubble.numberOfSuperbubbles() << "\n";
                     currentinsertTime = 0.0;
                 }
                 

                 

            }
            
            //timeFile << g.displayNumberEdges() << " " << g.displayNumberVertices() << "\n";
            
            if (mode == "perBatch") {
                start = std::chrono::high_resolution_clock::now();
                d.performDFS(&g);
                end = std::chrono::high_resolution_clock::now();
                elapsed = end - start;
                dfsTime = elapsed.count();
            
                start = std::chrono::high_resolution_clock::now();
                d.createMinI(&g); 
                end = std::chrono::high_resolution_clock::now();
                elapsed = end - start;
                minITime = elapsed.count();
            
                start = std::chrono::high_resolution_clock::now();
                d.createMaxO(&g);
                end = std::chrono::high_resolution_clock::now();
                elapsed = end - start;
                maxOTime = elapsed.count();
            

                start = std::chrono::high_resolution_clock::now();
                superbubble.computeSB(&d,&g);
                end = std::chrono::high_resolution_clock::now();
                elapsed = end - start;
                superbubbleTime = elapsed.count();
                totalTime = insertTime + dfsTime + minITime + maxOTime + superbubbleTime;
                
                timeFile << update_counter << " " << clean_hap_name << " " << insertTime << " " << dfsTime << " " << minITime << " " << maxOTime << " " << superbubbleTime << " " << totalTime << " " << superbubble.numberOfSuperbubbles() << "\n";

            }
            timeFile << g.displayNumberEdges() << " " << g.displayNumberVertices() << "\n";
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
        timeFile << "PeakMemory(KB): " << getMemoryUsageInKB() << "\n";
        timeFile.close();
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
    SimpleSuperbubbleH1 simple;
    simple.inputAndCompute(argv[1], argv[2]);
    return 0; 
}