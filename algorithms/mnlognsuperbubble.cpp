#include <sys/resource.h>
#include"superbubbleAlgorithms.h"
#include<fstream>
#include<chrono>
class IncrementalSuperbubbleDifferentAVLTree: public AlgorithmForSuperbubbles {
public:  
    long getMemoryUsageInKB() {
        struct rusage usage;
        getrusage(RUSAGE_SELF, &usage);
        return usage.ru_maxrss; // memory in kilobytes
    }   
    void inputAndCompute() {
        int m,x,y,z, cases;
        z=1;
        std::cin >> m;
        //std::ofstream outFile("/home/dell/Research/IncSuperbubble/output/Yeast/superbubbles/SuperbubbleincrementalDFSAVLTreeOutput.txt");
        std::ofstream timeFile("/home/user/Superbubble/IncrementalSuperbubble/output/Human/MiniGraph/timeFile/GRCh38-464.gfa/SuperbubbleincrementalDFSAVLTreeOutput.txt");
        AdjacencyListGraph adl(m) ;
        GraphContext g(&adl);
        IncrementalDFSTREE incdfs(&g);
        DFSCONTEXT d(&incdfs);
        SuperbubbleCalculation superbubble;
        
        if (timeFile.is_open()) {
              timeFile <<"             \n"; 
           } else {
               std::cerr << "Error opening time file at: " << std::endl;
        }
                
        for (int i = 1; i <= m; i++) {
           std::cin >> x >> y;
           std::cout << x <<" "<< y << std::endl;
           auto start = std::chrono::high_resolution_clock::now();  // Start time
           superbubble.insertE(&d,&g,x,y);
           superbubble.computeSB(&d,&g);
           
           auto end = std::chrono::high_resolution_clock::now();  // End time
           std::chrono::duration<double> elapsed = end - start;
           int superbubblecount = superbubble.numberOfSuperbubbles();
           if (timeFile.is_open()) {
              timeFile << z << " " << elapsed.count() << " " << superbubblecount << "\n" ;
           } else {
               std::cerr << "Error opening time file at: " << std::endl;
           }
           z++;
           
           
        } 

        timeFile << "PeakMemory(KB): " << getMemoryUsageInKB() << "\n";
        timeFile.close();
        std::fstream timeFileUpdate("/home/user/Superbubble/IncrementalSuperbubble/output/Human/MiniGraph/timeFile/GRCh38-464.gfa/SuperbubbleincrementalDFSAVLTreeOutput.txt", std::ios::in | std::ios::out);
        if (timeFileUpdate.is_open()) {
            timeFileUpdate.seekp(0, std::ios::beg);  // Move to the beginning of the file
            timeFileUpdate << g.displayNumberEdges() << " " << g.displayNumberVertices();
        } else {
            std::cerr << "Error opening time file at: " << std::endl;
        }
        
        /*  
        g.displayOutgoingGraph();
        g.displayIncomingGraph();
        std::cout << "Number of edges = " << g.displayNumberEdges() << std::endl;
        std::cout << "Number of vertices = " << g.displayNumberVertices() << std::endl;
        d.displayReverseTopologicalOrder();
        d.displayPostorder();
        d.displayReversePostorder();
        d.displayTopologicalOrder();
        d.displayMinI();
        d.displayMaxO();
        superbubble.printSuperbubbles();
        std::cout << "number of superbubbles are" << std::endl;
        std::cout << superbubble.numberOfSuperbubbles();
        */
        /*
        if (outFile.is_open()) {
            for (const pii& elem : superbubble.SB) {
                outFile << elem.first << " " << elem.second << "\n";
            }
            std::cout << "List written to file successfully at: "<< std::endl;
            outFile.close();
        } else {
            std::cerr << "Error opening file at: " << std::endl; 
        }
        
        */ 
    }
};
int main(){
	IncrementalSuperbubbleDifferentAVLTree incremental;
    incremental.inputAndCompute();
}
