#include"../lib/graph/contextGraph.h"
#include"../lib/DFS/contextdfs.h"
#include"../lib/graph/adjacencyListGraph.h"
#include"../lib/DFS/incrementalDFSAVLTree.h"
#include<fstream>
#include<string>

void incDFSAVLMiIMaxO(std::string& filename) { 
    int m,x,y; 
    
    std::string baseDirOutput = "/home/dell/Research/IncSuperbubble/output/HunderedRandomGraphs/IncrementalDFSAVLTree/";
    std::string baseDirInput = "/home/dell/Research/IncSuperbubble/inputGraphs/hunderedRandomGraphs/";
    std::string extension = ".txt";
    std::string outputFileName = baseDirOutput + filename + "incrementalDFSAvlTree.txt";
    filename = baseDirInput + filename + extension;
    std::ifstream infile(filename);
    std::ofstream outFile(outputFileName);
    infile >> m;
    AdjacencyListGraph adl(m) ;
    GraphContext g(&adl);
    IncrementalDFSTREE incdfs(&g);
    DFSCONTEXT d(&incdfs);
    for (int i = 1; i <= m; i++) {
           infile >> x >> y;
           std::cout << x << " "<<y<< std::endl;
           d.updateEdge(&g,x,y);
    }   
    
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
    
    
    if (outFile.is_open()) {
        // Write vector elements to file
        int i=0;
        outFile << "postorder" << std::endl;
        for (int element : d.dfsStrategy->postorder) {  
            outFile << i << " " << element << std::endl;
            i++;
        }
        i=0;
        outFile << "reversePostorder" << std::endl;
        for (int element : d.dfsStrategy->reversePostorder) {            
            outFile << i << " " << element << std::endl;
            i++;
        }
        i=0;
        outFile << "reverseTopologicalOrder" << std::endl;
        for (int element : d.dfsStrategy->reverseTopologicalOrder) {              
            outFile << i << " " << element << std::endl;
            i++;
        }
        i=0;
        outFile << "topOrd" << std::endl;
        for (int element : d.dfsStrategy->topOrd) {             
            outFile << i << " " << element << std::endl;
            i++;
        }
        i=0;
        outFile << "minI" << std::endl;
        for (int element : d.dfsStrategy->minI) {  
            outFile << i << " " << element << std::endl;
            i++;
        }
        i=0;
        outFile << "maxO" << std::endl; 
        for (int element : d.dfsStrategy->maxO) {
             outFile << i << " " << element << std::endl;
            i++;
        }
        outFile.close(); // Close file
        std::cout << "Vector written to file successfully!" << std::endl;
    } else {
        std::cerr << "Error opening file!" << std::endl;
    }
    
}



int main() {
    
    std::string filename;
    
    for(int i=1; i <= 100; i++) {
        filename = "randomGraph" + std::to_string(i);
        incDFSAVLMiIMaxO(filename);
    }
    
    //filename = "randomgraph2processed";
    //incDFSAVLMiIMaxO(filename);
	return 0;
}
