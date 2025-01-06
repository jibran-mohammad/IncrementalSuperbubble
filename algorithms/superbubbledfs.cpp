#include"superbubbleAlgorithms.h"
#include<fstream>
#include<chrono>
class SimpleSuperbubble: public AlgorithmForSuperbubbles {
public:    
    void inputAndCompute() {
        int m,x,y,z;
        z=1;
        std::cin >> m;
        //std::ofstream outFile("/home/dell/Research/IncSuperbubble/output/Yeast/superbubbles/SuperbubblesimpleDFSOutput.txt");
        std::ofstream timeFile("/home/dell/Research/IncSuperbubble/output/Yeast/timeFile/SuperbubbleSimpleTimeOutput.txt");
        AdjacencyListGraph adl(m) ;
        GraphContext g(&adl);
        SimpleDFS simdfs(&g);
        DFSCONTEXT d(&simdfs);
        SuperbubbleCalculation superbubble;
        
        
        if (timeFile.is_open()) {
              timeFile <<"             \n";
           } else {
               std::cerr << "Error opening time file at: " << std::endl;
        }
        
    
        d.updateEdge(&g,0,1);
        for (int i = 1; i <= m; i++) {
           std::cin >> x >> y;
           std::cout << x << " "<< y << std::endl;
           auto start = std::chrono::high_resolution_clock::now();  // Start time
           superbubble.insertE(&d,&g,x, y);
           d.performDFS(&g);
           d.createMinI(&g);
           d.createMaxO(&g);
           superbubble.computeSB(&d,&g);
           auto end = std::chrono::high_resolution_clock::now();  // End time
           std::chrono::duration<double> elapsed = end - start;
           
           if (timeFile.is_open()) {
              timeFile << z << " " << elapsed.count() << "\n" ;
           } else {
               std::cerr << "Error opening time file at: " << std::endl;
           }
           z++;
           
       }
       
       timeFile.close();
       
       std::fstream timeFileUpdate("/home/dell/Research/IncSuperbubble/output/Yeast/timeFile/SuperbubbleSimpleTimeOutput.txt", std::ios::in | std::ios::out);
       
       if (timeFileUpdate.is_open()) {
            timeFileUpdate.seekp(0, std::ios::beg);  // Move to the beginning of the file
            timeFileUpdate << g.displayNumberEdges() << " " << g.displayNumberVertices();
        } else {
            std::cerr << "Error opening time file at: " << std::endl;
        }
        
       /*
       g.displayOutgoingGraph();
       g.displayIncomingGraph();
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
	SimpleSuperbubble simple;
    simple.inputAndCompute();
}