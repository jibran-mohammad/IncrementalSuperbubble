#include "baseDFS.h"
    
    void BaseDFS::dfs(GraphContext *G) { 
        //this->a=0;
        //this->b=0;
        #ifndef BENCHMARK
        auto start = std::chrono::high_resolution_clock::now();
        auto end = std::chrono::high_resolution_clock::now();   
        std::chrono::duration<double> elapsed = end - start;
        start = std::chrono::high_resolution_clock::now();
        #endif
        for (int i = 0; i <this->n; ++i) {
            if (!this->visit[i]) {
               this->dfsUtil(G, i);
           }
        }
        #ifndef BENCHMARK
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        this->latestTimes.dfsTime = elapsed.count();
        this->latestTimes.dfsTime = this->latestTimes.dfsTime - this->latestTimes.topologicalOrderTime;
        #endif

       //std::cout <<" Nodes and edges " << this->a << " " << this->b << std::endl;
    }
 
    void BaseDFS::dfsUtil(GraphContext *G, int v) {
        #ifndef BENCHMARK
        auto start = std::chrono::high_resolution_clock::now();
        auto end = std::chrono::high_resolution_clock::now();   
        std::chrono::duration<double> elapsed = end - start;
        #endif
        this->visit[v] = 1; // Mark the current vertex as visited
        //this->a++;
        // Recur for all adjacent vertices of the current vertex
        for (int u : G->displayOutgoingNeighbours(v)) {
            //this->b++;
            if (!this->visit[u]) {
                this->dfsUtil(G, u);
            }
        } 
          
        // After visiting all adjacent vertices, add the current vertex to the postorder
        #ifndef BENCHMARK
        start = std::chrono::high_resolution_clock::now();
        #endif
        this->postorder[v]=this->postorderNumber;
        this->reversePostorder[this->postorderNumber] = v;
        this->reverseTopologicalOrder[v]=(this->n-1) - this->postorder[v];
        this->topOrd[this->reverseTopologicalOrder[v]]=v;
        this->postorderNumber++; 
        #ifndef BENCHMARK
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        this->latestTimes.currenttopologicalOrderTime = elapsed.count();
        this->latestTimes.topologicalOrderTime += elapsed.count();  
        #endif
        
         
    }

    void BaseDFS::printTopologicalOrder() const {
    	std::cout << "The topological order of the given graph is " <<std::endl;
    	for(int i=0; i<this->n; i++) {
    		std::cout << i << "->" << this->topOrd[i];
    		std::cout << std::endl;
    	}

    }

    void BaseDFS::printreversePostOrder() const {
        std::cout << "The ReversePostorder order of the given graph is " <<std::endl;
        for(int i=0; i<this->n; i++) {
            std::cout << i << "->" << this->reversePostorder[i];
            std::cout << std::endl;
        }

    }

    void BaseDFS::printPostOrder() const {
        std::cout << "The Postorder order of the given graph is " <<std::endl;
        
        for(int i=0; i<this->n; i++) {
            std::cout << i << "->" << this->postorder[i];
            std::cout << std::endl;
        }
    }

    void BaseDFS::printReverseTopologicalOrder() const {
        std::cout << "The reverse topological order of the given graph is " <<std::endl;
        for(int i=0; i<this->n; i++) {
            std::cout << i << "->" << this->reverseTopologicalOrder[i];
            std::cout << std::endl;
        }
    }
    /*
      this function is used for inserting the edge. It figures out the number of whether any new nodes are
      added or not, it added, then it resizes all the datastructures accordingly and calls reset, otherwise
      it just calls reset and makes postordernumber 0.
    */
    void BaseDFS::insert(GraphContext *g, int u, int v) {
        int x;
        x=this->n;
        g->insertEdge(u,v);
        this->n=g->displayNumberVertices();
        x=this->n-x;
        if(x>0) {
            this->postorderNumber = 0;
            this->reversePostorder.resize(this->n);
            this->oldReverseTopologicalOrder.resize(this->n);
            this->postorder.resize(this->n);
            this->visit.resize(this->n);
            this->topOrd.resize(this->n);
            this->reverseTopologicalOrder.resize(this->n);
            this->minI.resize(this->n);
            this->maxO.resize(this->n);
            this->reset();
        }
        else {
            this->postorderNumber = 0;
            this->reset();
        }
    }

    void BaseDFS::reset() {
        std::fill(this->reversePostorder.begin(), this->reversePostorder.end(), 0);
        std::fill(this->postorder.begin(), this->postorder.end(), 0);
        std::fill(this->visit.begin(), this->visit.end(), 0);
        std::fill(this->topOrd.begin(), this->topOrd.end(), 0);
        std::fill(this->reverseTopologicalOrder.begin(), this->reverseTopologicalOrder.end(), 0);
        std::fill(this->oldReverseTopologicalOrder.begin(), this->oldReverseTopologicalOrder.end(), 0);
        std::fill(this->minI.begin(), this->minI.end(), 0);
        std::fill(this->maxO.begin(), this->maxO.end(), 0);
    }
