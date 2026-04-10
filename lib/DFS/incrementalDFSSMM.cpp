#include"incrementalDFSSMM.h"

    
    /*
      insert function used for insertion of an edge and for maintaining all relevant datastructures like 
      reversetopological order, minI, maxO etc. It handles five cases. 
      1) edge (u,v) inserted, one new node inserted and u is new node.
      2) edge (u,v) inserted, one new node inserted and v is new node.
      3) edge (u,v) inserted, two new nodes inserted.
      4) edge (u,v) inserted, edge is from right to left.
      5) edge (u,v) inserted, edge if from left to right.  
    */
    void IncrementalDFSSMM::insert(GraphContext *g, int u, int v) {
        int x;
        #ifndef BENCHMARK
        this->latestTimes.currentinsertTime = 0.0;
        this->latestTimes.currentdfsTime = 0.0;
        this->latestTimes.currenttopologicalOrderTime = 0.0;
        this->latestTimes.currentminITime = 0.0;
        this->latestTimes.currentmaxOTime = 0.0;
        this->edgeTypes.currentcase1 = 0;
        this->edgeTypes.currentcase2 = 0;
        this->edgeTypes.currentcase3 = 0;
        this->edgeTypes.currentcase4 = 0;
        this->edgeTypes.currentcase5 = 0;
        this->vertexCount.currentcountAffectedVertices = 0;
        this->vertexCount.currentcountExploredVertices = 0;
        this->vertexCount.currentTotalEdges = 0;

        auto start = std::chrono::high_resolution_clock::now();
        auto end = std::chrono::high_resolution_clock::now();   
        std::chrono::duration<double> elapsed = end - start;

        start = std::chrono::high_resolution_clock::now();
        #endif
        x=this->n;
        g->insertEdge(u,v);
        this->n=g->displayNumberVertices();
        this->vertexCount.currentTotalEdges = g->displayNumberEdges();
        this->vertexCount.totalEdges += this->vertexCount.currentTotalEdges;
        x=this->n-x;
        //std::cout << "original U= "<< u << " " << "original v= " << v << std::endl;
        this->u=g->strategy->reverseMapping[u];
        this->v=g->strategy->reverseMapping[v];
        u=g->strategy->reverseMapping[u];
        v=g->strategy->reverseMapping[v];
        #ifndef BENCHMARK
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        this->latestTimes.currentinsertTime = elapsed.count();
        this->latestTimes.insertTime += elapsed.count();
        #endif
        //std::cout << "new U= "<< u << " " << "new v= " << v << std::endl;
        if(x>0) {
            //std::cout << " x= " << x << std::endl;
            //std::cout << " n= " << n << std::endl;
            #ifndef BENCHMARK
            start = std::chrono::high_resolution_clock::now();
            #endif
            this->reversePostorder.resize(this->n);
            this->postorder.resize(this->n);
            this->visit.resize(this->n);
            this->topOrd.resize(this->n);
            this->reverseTopologicalOrder.resize(this->n);
            this->minI.resize(this->n);
            this->maxO.resize(this->n);
            this->marked.resize(this->n);
            this->intermediatepostorder.resize(this->n);
            this->postorderChanged.resize(this->n);
            #ifndef BENCHMARK
            end = std::chrono::high_resolution_clock::now();
            elapsed = end - start;
            this->latestTimes.currentinsertTime += elapsed.count();
            this->latestTimes.insertTime += elapsed.count();
            #endif
            //std::cout << "after resize" << std::endl;
            if(x == 1) { 
               //std::cout << "x is equal to 1" << std::endl;
               //case where u is new node
               if(g->strategy->newu == true){
                   this->cases = 1;
                   #ifndef BENCHMARK
                   this->edgeTypes.case1 +=1;
                   this->edgeTypes.currentcase1 = 1;
                   #endif
                   //std::cout << "case 1" << std::endl;

                   #ifndef BENCHMARK
                   start = std::chrono::high_resolution_clock::now();
                   #endif
                   this->postorder[u] = this->postorder[0];
                   this->reversePostorder[postorder[u]] = u;
                   this->reverseTopologicalOrder[u] = (this->n - 1) -  this->postorder[u];
                   this->topOrd[reverseTopologicalOrder[u]] = u;           
                   
                   this->postorder[0] = this->postorder[u] + 1;
                   this->reversePostorder[postorder[0]] = 0;
                   this->reverseTopologicalOrder[0] = (this->n - 1) -  this->postorder[0];
                   this->topOrd[reverseTopologicalOrder[0]] = 0;                     

                    for(int i=0; i<this->n; i++) {
                        this->reverseTopologicalOrder[i] = (this->n - 1) - postorder[i];
                        this->topOrd[reverseTopologicalOrder[i]] = i;
                    }
                    #ifndef BENCHMARK
                    end = std::chrono::high_resolution_clock::now();
                    elapsed = end - start;
                    this->latestTimes.currenttopologicalOrderTime += elapsed.count();
                    this->latestTimes.topologicalOrderTime += elapsed.count();
                    #endif
                    g->strategy->newu = false;
                }    
                else {
                    //std::cout << "case 2" << std::endl;
                    this->cases=2;
                    #ifndef BENCHMARK
                    this->edgeTypes.case2 +=1;
                    this->edgeTypes.currentcase2 = 1;


                    start = std::chrono::high_resolution_clock::now();
                    #endif
                    this->postorder[v] = this->postorder[u];
                    for(int i=postorder[u]; i < this->n - 1; i++) {
                      this->postorder[reversePostorder[i]]++;
                    }
                    
                    for(int i=0;i < this->n;i++) {
                      this->reversePostorder[this->postorder[i]]=i;
                      this->reverseTopologicalOrder[i]=(this->n) - 1 - this->postorder[i];
                      this->topOrd[this->reverseTopologicalOrder[i]] =i; 
                    }
                    #ifndef BENCHMARK
                    end = std::chrono::high_resolution_clock::now();
                    elapsed = end - start;
                    this->latestTimes.currenttopologicalOrderTime += elapsed.count();
                    this->latestTimes.topologicalOrderTime += elapsed.count();
                    #endif

                    g->strategy->newv= false;    
                }    
            }
            else {
               this->cases = 3;
               #ifndef BENCHMARK
               this->edgeTypes.case3 +=1;
               this->edgeTypes.currentcase3 = 1;
               start = std::chrono::high_resolution_clock::now();
               #endif
               //std::cout << "if x is 2" << std::endl; 
               this->postorder[v] = this->postorder[0];
               this->reversePostorder[postorder[v]] = v;
               this->reverseTopologicalOrder[v] = (this->n - 1) -  this->postorder[v];
               this->topOrd[reverseTopologicalOrder[v]] = v;
            
               this->postorder[u] = this->postorder[v] + 1;
               this->reversePostorder[postorder[u]] = u;
               this->reverseTopologicalOrder[u] = (this->n - 1) -  this->postorder[u];
               this->topOrd[reverseTopologicalOrder[u]] = u;            

               this->postorder[0] = this->postorder[u] + 1;
               this->reversePostorder[postorder[0]] = 0;
               this->reverseTopologicalOrder[0] = (this->n - 1) -  this->postorder[0];
               this->topOrd[reverseTopologicalOrder[0]] = 0;
               for(int i=0; i<this->n;i++) {
                   this->reverseTopologicalOrder[i] = (this->n - 1) - postorder[i];
                   this->topOrd[reverseTopologicalOrder[i]] = i;
                }
                #ifndef BENCHMARK
                end = std::chrono::high_resolution_clock::now();
                elapsed = end - start;
                this->latestTimes.currenttopologicalOrderTime += elapsed.count();
                this->latestTimes.topologicalOrderTime += elapsed.count();
                #endif

               g->strategy->newu = false;
               g->strategy->newv = false;
                           
            }
        }    
        else {
            
            //std::cout << "if no new node added" << std::endl;
            if(this->postorder[u] < this->postorder[v]) {
                this->cases = 4;
                this->edgeTypes.case4 +=1;
                this->edgeTypes.currentcase4 = 1;
                #ifndef BENCHMARK
                start = std::chrono::high_resolution_clock::now();
                #endif
                this->partialDFS(g,u,v);
                #ifndef BENCHMARK
                end = std::chrono::high_resolution_clock::now();
                elapsed = end - start;
                this->latestTimes.currentdfsTime += elapsed.count();
                this->latestTimes.dfsTime += elapsed.count();
                #endif
    		    this->update(g);
                #ifndef BENCHMARK
                start = std::chrono::high_resolution_clock::now();
                #endif
                this->reset();
                #ifndef BENCHMARK
                end = std::chrono::high_resolution_clock::now();
                elapsed = end - start;
                this->latestTimes.currenttopologicalOrderTime += elapsed.count();
                this->latestTimes.topologicalOrderTime += elapsed.count();
                #endif

    	    }
            else {
                this->cases = 5;
                #ifndef BENCHMARK
                this->edgeTypes.case5 +=1;
                this->edgeTypes.currentcase5 = 1;
                #endif
            }   
        }
        /*
        g->displayOutgoingGraph();
        g->displayIncomingGraph();
        this->printPostOrder();
        this->printreversePostOrder();
        this->printReverseTopologicalOrder();
        this->printTopologicalOrder();
        std::cout << "incoming tree is " << std::endl;
        this->printtreeminIncoming();
        std::cout << "outgoing tree is " << std::endl;
        this->printtreemaxOutgoing();
        */
   }

    void IncrementalDFSSMM::partialDFS(GraphContext *g,int w, int z) {
    	if(this->marked[z] == 0 && (this->postorder[z] > this->postorder[this->u])) {
            this->marked[z] =1;
            for (int u : g->displayOutgoingNeighbours(z)) {
                this->partialDFS(g,z,u);
            }
            this->m = this->m+1;
            this->postorderChanged[z] = m;
     
    	}
    }

    void IncrementalDFSSMM::update(GraphContext *g) {
        int shift = this->m;
    	int index;
        
        #ifndef BENCHMARK
        auto start = std::chrono::high_resolution_clock::now();
        auto end = std::chrono::high_resolution_clock::now();   
        std::chrono::duration<double> elapsed = end - start;


        start = std::chrono::high_resolution_clock::now();
        #endif
    	
        for(int i = this->postorder[this->u]; i <= this->postorder[this->v]; i++) {
            this->vertexCount.currentcountAffectedVertices ++;
            index = this->reversePostorder[i];
            if(this->marked[index]==1) {
            	this->intermediatepostorder[index] = this->postorder[this->u] + this->postorderChanged[index] - 1;
            	this->marked[index] = 0;
            }
            else {
            	this->intermediatepostorder[index] = this->postorder[this->u] + shift;
            	shift = shift + 1;
            }
    	}
        #ifndef BENCHMARK
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        this->latestTimes.currenttopologicalOrderTime += elapsed.count();
        this->latestTimes.topologicalOrderTime += elapsed.count();
        

        start = std::chrono::high_resolution_clock::now(); 
        #endif
        for(int i=0; i<this->n; i++) {
            if(this->intermediatepostorder[i]!=0) {
                //std::cout << "i = " << i << std::endl;
                if(i == 0 || i == 1) {
                    continue;
                }       
            	this->postorder[i] = this->intermediatepostorder[i];
            	this->reversePostorder[this->intermediatepostorder[i]]=i;
                this->reverseTopologicalOrder[i]=(this->n-1) - this->postorder[i];
                this->topOrd[this->reverseTopologicalOrder[i]] = i;
            }    
        }
        #ifndef BENCHMARK
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        this->latestTimes.currenttopologicalOrderTime += elapsed.count();
        this->latestTimes.topologicalOrderTime += elapsed.count();
        this->vertexCount.countAffectedVertices += this->vertexCount.currentcountAffectedVertices;
        #endif
    }
    
    void IncrementalDFSSMM::reset() {
        std::fill(this->marked.begin(), this->marked.end(), 0);
        std::fill(this->postorderChanged.begin(), this->postorderChanged.end(), 0);
        std::fill(this->intermediatepostorder.begin(), this->intermediatepostorder.end(), 0);
        this->m=0;
    }

    void IncrementalDFSSMM::calculateMaxO(GraphContext *g) {
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

    void IncrementalDFSSMM::calculateMinI(GraphContext *g) {
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
    void IncrementalDFSSMM::printMinI() const {
        std::cout << "MinI of the graph is" << std::endl;
        for(int i=0; i < this->n; i++) {
            std::cout << i << " -> " << this->minI[i] << std::endl; 
        }
    }

    void IncrementalDFSSMM::printMaxO() const {
        std::cout << "MaxO of the graph is " << std::endl;
        for(int i=0; i< this->n; i++) {
            std::cout << i << " -> " << this->maxO[i] << std::endl;
        }
    }
