#include"simpleDFS.h"
    /*
       this function calculate the minimum incoming edge of every node, it will take O(m+n) time, because it
       is exploring every edge of the graph, except for node 0 and 1, node 0 is the artifical source and its 
       minimum incoming edge will be infinity, and minimum incoming edge of 1 would be 0, because edge 0,1
       is artificially inserted. Datastructure used for minI is simple array.
    */
    void SimpleDFS::calculateMinI(GraphContext *g) {
        int min = INT_MAX;
        int element = INT_MAX;
        this->minI[0] = INT_MAX;
        this->minI[1] = 0;
        for(int i=2; i < this->n; i++) {
            if(g->strategy->incoming[i] == 2) {
                element = 0;
            }
            else {
               for(int u: g->displayIncomingNeighbours(i)) {
                   if(u == 0){
                     continue;
                   }
                   if(this->reverseTopologicalOrder[u] < min ) {
                      min = this->reverseTopologicalOrder[u];
                      element = this->topOrd[min];
                   }
               }
            }   
            this->minI[i] = element;
            min = INT_MAX;
        }
    }
    
    /*
       this function calculate the maximum outgoin edge of every node, it will take O(m+n) time, because it
       is exploring every edge of the graph, except for node 0 and 1, node 0 is the artifical source and its 
       maximum outgoing edge will be 1,because edge 0,1 is artificially inserted. and maximum outgoing edge 
       of 1 would be infinity, because it doesn't outgoing edge. Datastructure used for maxO is simple array. 
    */
    void SimpleDFS::calculateMaxO(GraphContext *g) {
        int max = 0;
        this->maxO[1]= INT_MAX;
        this->maxO[0]= 1;
        int element = INT_MAX;
        for(int i =2; i < this->n; i++) {
            if(g->strategy->outgoing[i]==2) {
                element = 1;
            }
            else {
                for(int u: g->displayOutgoingNeighbours(i)) {
                    if(u == 1) {
                        continue;
                    }
                    if(this->reverseTopologicalOrder[u] > max) {
                        max = this->reverseTopologicalOrder[u];
                        element = this->topOrd[max];
                    }
                }
            }  
            this->maxO[i] = element;
            max = 0;
        }
    }

    void SimpleDFS::printMinI() const {
        std::cout << "MinI of the graph is" << std::endl;
        for(int i=0; i < this->n; i++) {
            std::cout << i << " -> " << this->minI[i] << std::endl; 
        }
    }

    void SimpleDFS::printMaxO() const {
        std::cout << "MaxO of the graph is " << std::endl;
        for(int i=0; i< this->n; i++) {
            std::cout << i << " -> " << this->maxO[i] << std::endl;
        }
    }
