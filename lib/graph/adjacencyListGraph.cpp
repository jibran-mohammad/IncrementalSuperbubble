#include <iostream>
#include "adjacencyListGraph.h"


// Define functions for the AdjacencyListGraph class


   
   
   /* 
   add's a single edge to the list
   */

   void AdjacencyListGraph::addEdge(int x, int y) {
       this->artificialParentRemoved = false;
       if (x < 0 || y < 0) {
           std::cout << "ERROR: Incorrect arguments addEdge " << x << " " << y << std::endl;
           return;
       }
       auto it = this->reverseMapping.find(x);
       if (it == reverseMapping.end()) {
           this->addNode(x);
           this->newu=true;
       }

       auto its = this->reverseMapping.find(y);
       
       if (its == reverseMapping.end()) {
           this->addNode(y);
           this->newv=true;           
       }

       //std::cout << x << " " << y << std::endl;
       //std::cout << this->reverseMapping[x] << " " << this->reverseMapping[y] << std::endl;

       if(this->outgoing[this->reverseMapping[x]] == 2) {
           this->outgoing[this->reverseMapping[x]] = 1;
       }
       else {
           this->outgoing[this->reverseMapping[x]] = 1;
       }

       auto t = std::find(this->adjLOutgoing[this->reverseMapping[x]].begin(),this->adjLOutgoing[this->reverseMapping[x]].end(),this->reverseMapping[y]);
       if(t == this->adjLOutgoing[this->reverseMapping[x]].end()) {
           this->adjLOutgoing[this->reverseMapping[x]].push_back(this->reverseMapping[y]);  // Add to outgoing edges
           this->adjLIncoming[this->reverseMapping[y]].push_back(this->reverseMapping[x]);  // Add to incoming edges (for completeness) 
       }
       if(this->incoming[this->reverseMapping[y]] == 2) {
           this->incoming[this->reverseMapping[y]] = 1 ;
           /*
              y had only the artificial parent 0 and now has a real one: delete
              the artificial edge 0->y. Otherwise the DFS from 0 keeps y's
              subtree under 0 instead of moving it under x, and the
              superbubble computation, which needs a DFS order, misses
              superbubbles. Artificial edges to the sink 1 are kept: they do
              not change the DFS tree (1 is the first vertex visited) and are
              ignored by maxO.
           */
           this->adjLOutgoing[0].remove(this->reverseMapping[y]);
           this->adjLIncoming[this->reverseMapping[y]].remove(0);
           this->sizeE--;
           this->artificialParentRemoved = true;
       }
       else {
           this->incoming[this->reverseMapping[y]] = 1;
       }
       /*
          Artificial edge added from 0 to x, if x is the new node(indicated by value at index(given by x) is zero).
          After adding artificial edge, the value at index(given by x) is made to be 2. 
       */
       if(this->incoming[this->reverseMapping[x]] == 0 && this->reverseMapping[x]!=0) {
          this->adjLOutgoing[0].push_back(this->reverseMapping[x]);
          this->sizeE++;
          this->adjLIncoming[this->reverseMapping[x]].push_back(0);
          this->incoming[this->reverseMapping[x]] = 2 ;
       }
       /*
          Artificial edge added from y to 1, if y is the new node(indicated by value at index(given by y) is zero).
          After adding artificial edge, the value at index(given by y) is made to be 2. 
       */
       if((this->outgoing[this->reverseMapping[y]] == 0) && (this->reverseMapping[y]!=1)) {
          this->adjLIncoming[1].push_back(this->reverseMapping[y]);
          this->sizeE++;
          this->adjLOutgoing[this->reverseMapping[y]].push_back(1);
          this->outgoing[this->reverseMapping[y]] = 2;
       }
       
          
    }
   
   
   /*
     Add's a single node to the graph.
   */
   void AdjacencyListGraph::addNode(int x) {
        
        this->reverseMapping[x] = sizeV;
        this->mapping.push_back(x);
        this->adjLOutgoing.push_back(std::list<int>{});
        this->adjLIncoming.push_back(std::list<int>{});
        sizeV++;
        this->incoming.push_back(0);
        this->outgoing.push_back(0);

        
   }


   /*
   printing the graph with  outgoing edges
   */

   void AdjacencyListGraph::printGOutgoing() const {
       std::cout << "Outgoing Graph" << std::endl;
       std::cout << this->sizeV << " " << this->sizeE << std::endl;

       for (int i = 0; i < this->sizeV; i++) {
           std::cout << i << ": ";
           for (int neighbor : this->adjLOutgoing[i]) {
               std::cout << neighbor << " ";
           }
           std::cout << std::endl;
       }
   }

   /*
   Printing the graph with incoming edges
   */

   void AdjacencyListGraph::printGIncoming() const {
       std::cout << "Incoming Graph" << std::endl;
       std::cout << this->sizeV << " " << this->sizeE << std::endl;

       for (int i = 0; i < this->sizeV; i++) {
           std::cout << i << ": ";
           for (int neighbor : this->adjLIncoming[i]) {
               std::cout << neighbor << " ";
           }
           std::cout << std::endl;
       }
   }
   
   /*
   getting the neighbours of the outgoing adjacency list
   */

   const std::list<int>& AdjacencyListGraph::getNeighboursOutgoing (int x) const {
       // Return a reference to the outgoing neighbors for consistency
       static const std::list<int> empty_list;  

       if (x < 0 || x >= (int)this->adjLOutgoing.size()) {  
           return empty_list;
       }
        return this->adjLOutgoing[x];
  
   }
   
   /*
   getting the neigbours of the incoming adjacency list
   */

   const std::list<int>& AdjacencyListGraph::getNeighboursIncoming(int x) const {
       static const std::list<int> empty_list;  

       if (x < 0 || x >= (int)this->adjLIncoming.size()) {  
           return empty_list;
       }
       return this->adjLIncoming[x];
       
   }     

   /*
   Returns the number of vertices in the graph
   */

   int AdjacencyListGraph::getNumberOfVertices() const {
       return this->sizeV;
   }

   /*
   Returns the number of edges in the graph
   */

   int AdjacencyListGraph::getNumberOfEdges() const {
       return this->sizeE;
   }

    int AdjacencyListGraph::getInDegree(int x) const {
        if (x < 0 || x >= (int)this->adjLIncoming.size()) {
            std::cout << "ERROR: Node " << x << " not found" << std::endl;
            return -1;
        }
        return this->adjLIncoming[x].size();
    }

    int AdjacencyListGraph::getOutDegree(int x) const {
        
        if (x < 0 || x >= (int)this->adjLOutgoing.size()) {
            std::cout << "ERROR: Node " << x << " not found" << std::endl;
            return -1;
        }
        return this->adjLOutgoing[x].size();
    }

    int AdjacencyListGraph::getDegree(int x) const {
        int in = getInDegree(x);
        int out = getOutDegree(x);
        if (in == -1 || out == -1) return -1;
        return in + out;
    }
