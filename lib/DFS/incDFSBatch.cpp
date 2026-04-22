#include"incDFSBatch.h"
    /*
      insert function used for insertion of an edge and for maintaining all relevant datastructures like 
      reversetopological order, minI, maxO etc. It handles five cases. 
      1) edge (u,v) inserted, one new node inserted and u is new node.
      2) edge (u,v) inserted, one new node inserted and v is new node.
      3) edge (u,v) inserted, two new nodes inserted.
      4) edge (u,v) inserted, edge is from right to left.
      5) edge (u,v) inserted, edge if from left to right.
      the difference from the O(m^2logn) is the number of vertices we need to explore whose topological order
      got changed due to insertion of an edge. We are only changing the vertices in the special balanced binary
      search tree whos reachability changed.  
    */

    void IncrementalDFSTREEBatch::insert(GraphContext *g, int u, int v) {
        

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
            this->vertexCount.currentTotalAffectedEdges = 0;
            this->vertexCount.currentTotalExploredEdges = 0;

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
                this->reachable.resize(this->n);
                this->globalReachable.resize(this->n);
                this->minI.resize(this->n);
                this->maxO.resize(this->n);
                this->minIbst.resize(this->n);
                this->maxObst.resize(this->n);
                this->marked.resize(this->n);
                this->intermediatepostorder.resize(this->n);
                this->postorderChanged.resize(this->n);
                #ifndef BENCHMARK
                end = std::chrono::high_resolution_clock::now();
                elapsed = end - start;
                this->latestTimes.currentinsertTime += elapsed.count();
                this->latestTimes.insertTime += elapsed.count();
                #endif

                if(x == 1) {
                   //std::cout << "x is equal to 1" << std::endl;
                   //case where u is new node
                    this->cases = 1;
                   if(g->strategy->newu == true){
                       //std::cout << "case 1" << std::endl;
                       #ifndef BENCHMARK 
                       this->edgeTypes.case1 +=1;
                       this->edgeTypes.currentcase1 =1;

                       //std::cout << "case 1" << std::endl;
                       start = std::chrono::high_resolution_clock::now();
                       #endif
                       DifferentBalancedBinaryTree binaryTreeIncoming(this->reverseTopologicalOrder);
                       end = std::chrono::high_resolution_clock::now();
                       #ifndef BENCHMARK
                       elapsed = end - start;
                       this->latestTimes.currentminITime += elapsed.count();
                       this->latestTimes.minITime += elapsed.count();

                       start = std::chrono::high_resolution_clock::now();
                       #endif
                       DifferentBalancedBinaryTree binaryTreeOutgoing(this->reverseTopologicalOrder);
                       #ifndef BENCHMARK
                       end = std::chrono::high_resolution_clock::now();
                       elapsed = end - start;
                       this->latestTimes.currentmaxOTime += elapsed.count();
                       this->latestTimes.maxOTime += elapsed.count();

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

                       
                    
                       for(int i=0; i<this->n;i++) {
                          this->reverseTopologicalOrder[i] = (this->n - 1) - postorder[i];
                          this->topOrd[reverseTopologicalOrder[i]] = i;
                       }
                       #ifndef BENCHMARK
                       end = std::chrono::high_resolution_clock::now();
                       elapsed = end - start;
                       this->latestTimes.currenttopologicalOrderTime += elapsed.count();
                       this->latestTimes.topologicalOrderTime += elapsed.count();

                       start = std::chrono::high_resolution_clock::now();
                       #endif
                       this->maxObst[u]= binaryTreeOutgoing;
                       this->maxObst[u].addNode(v,0);
                       this->maxO[u] =  this->maxObst[u].findMax();
                       #ifndef BENCHMARK
                       end = std::chrono::high_resolution_clock::now();
                       elapsed = end - start;
                       this->latestTimes.currentmaxOTime += elapsed.count();
                       this->latestTimes.maxOTime += elapsed.count();

                       start = std::chrono::high_resolution_clock::now();
                       #endif
                       this->minIbst[u]= binaryTreeIncoming;
                       this->minIbst[v].addNode(u,0);
                       this->minI[v] = this->minIbst[v].findMin();
                       #ifndef BENCHMARK
                       end = std::chrono::high_resolution_clock::now();
                       elapsed = end - start;
                       this->latestTimes.currentminITime += elapsed.count();
                       this->latestTimes.minITime+= elapsed.count();
                       #endif
                       g->strategy->newu = false;
                   }    
                   else {
                        //v is the new node
                        //std::cout << "case 2" << std::endl;
                        this->cases=2;
                        #ifndef BENCHMARK
                        this->edgeTypes.case2 +=1;
                        this->edgeTypes.currentcase2 =1;
                        start = std::chrono::high_resolution_clock::now();
                        #endif
                        DifferentBalancedBinaryTree binaryTreeIncoming(this->reverseTopologicalOrder);
                        #ifndef BENCHMARK
                        end = std::chrono::high_resolution_clock::now();
                        elapsed = end - start;
                        this->latestTimes.currentminITime += elapsed.count();
                        this->latestTimes.minITime += elapsed.count();

                        start = std::chrono::high_resolution_clock::now();
                        #endif
                        DifferentBalancedBinaryTree binaryTreeOutgoing(this->reverseTopologicalOrder);
                        #ifndef BENCHMARK
                        end = std::chrono::high_resolution_clock::now();
                        elapsed = end - start;
                        this->latestTimes.currentmaxOTime += elapsed.count();
                        this->latestTimes.maxOTime += elapsed.count();

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

                        start = std::chrono::high_resolution_clock::now();
                        #endif
                        binaryTreeIncoming.addNode(u,0);
                        this->minIbst[v] = binaryTreeIncoming;
                        this->minI[v] = this->minIbst[v].findMin();
                        #ifndef BENCHMARK
                        end = std::chrono::high_resolution_clock::now();
                        elapsed = end - start;
                        this->latestTimes.currentminITime += elapsed.count();
                        this->latestTimes.minITime += elapsed.count();
                        
                        start = std::chrono::high_resolution_clock::now();
                        #endif
                        this->maxObst[u].addNode(v,0);
                        this->maxObst[v] = binaryTreeOutgoing;
                        this->maxO[u] = this->maxObst[u].findMax();
                        this->maxO[v] = 1;
                        #ifndef BENCHMARK
                        end = std::chrono::high_resolution_clock::now();
                        elapsed = end - start;
                        this->latestTimes.currentmaxOTime += elapsed.count();
                        this->latestTimes.maxOTime += elapsed.count();
                        #endif
                        g->strategy->newv= false;
                    }
                }
                else {
                   //std::cout << "if x is 2" << std::endl; 
                   this->cases = 3;
                   #ifndef BENCHMARK
                   this->edgeTypes.case3 +=1;
                   this->edgeTypes.currentcase3 =1;
                   start = std::chrono::high_resolution_clock::now();
                   #endif
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

                   start = std::chrono::high_resolution_clock::now();
                   #endif
                   DifferentBalancedBinaryTree binaryTreeIncomingu(this->reverseTopologicalOrder);
                   DifferentBalancedBinaryTree binaryTreeIncomingv(this->reverseTopologicalOrder);
                   binaryTreeIncomingv.addNode(u,0);
                   this->minIbst[u]= binaryTreeIncomingu;
                   this->minIbst[v]= binaryTreeIncomingv;
                   this->minI[u] = 0;
                   this->minI[v] = u;
                   #ifndef BENCHMARK
                   end = std::chrono::high_resolution_clock::now();
                   elapsed = end - start;
                   this->latestTimes.currentminITime += elapsed.count();
                   this->latestTimes.minITime += elapsed.count();

                   start = std::chrono::high_resolution_clock::now();
                   #endif
                   DifferentBalancedBinaryTree binaryTreeOutgoingu(this->reverseTopologicalOrder);
                   DifferentBalancedBinaryTree binaryTreeOutgoingv(this->reverseTopologicalOrder);
                   binaryTreeOutgoingu.addNode(v,0);
                   this->maxObst[u]= binaryTreeOutgoingu;
                   this->maxObst[v]= binaryTreeOutgoingv;
                   this->maxO[u] = v;
                   this->maxO[v] = 1;
                   #ifndef BENCHMARK
                   end = std::chrono::high_resolution_clock::now();
                   elapsed = end - start;
                   this->latestTimes.currentmaxOTime += elapsed.count();
                   this->latestTimes.maxOTime += elapsed.count(); 
                   #endif
                   g->strategy->newu = false;
                   g->strategy->newv = false;
                   }           
                }
            else {
                
                //std::cout << "if no new node added" << std::endl;
                if(this->postorder[u] < this->postorder[v]) {
                    //std::cout << "left to right edge" << std::endl;
                    this->cases = 4;
                    #ifndef BENCHMARK
                    this->edgeTypes.case4 +=1;
                    this->edgeTypes.currentcase4 =1;
                    start = std::chrono::high_resolution_clock::now();
                    #endif
                    this->partialDFS(g,u,v);
                    #ifndef BENCHMARK
                    end = std::chrono::high_resolution_clock::now();
                    elapsed = end - start;
                    this->latestTimes.currentdfsTime += elapsed.count();
                    this->latestTimes.dfsTime += elapsed.count();
                    this->vertexCount.totalExploredEdges += this->vertexCount.currentTotalExploredEdges;
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
                    //std::cout << "right to left  edge" << std::endl;
                    this->cases = 5;
                    #ifndef BENCHMARK
                    this->edgeTypes.case5 +=1;
                    this->edgeTypes.currentcase5 =1;
                    start = std::chrono::high_resolution_clock::now();
                    #endif
                    this->maxObst[u].addNode(v,0);
                    this->maxO[u] = (this->maxObst[u].findMax());
                    #ifndef BENCHMARK
                    end = std::chrono::high_resolution_clock::now();
                    elapsed = end - start;
                    this->latestTimes.currentmaxOTime += elapsed.count();
                    this->latestTimes.maxOTime += elapsed.count();
     
                    start = std::chrono::high_resolution_clock::now();
                    #endif
                    this->minIbst[v].addNode(u,0);
                    this->minI[v] = (this->minIbst[v].findMin());
                    #ifndef BENCHMARK
                    end = std::chrono::high_resolution_clock::now();
                    elapsed = end - start;
                    this->latestTimes.currentminITime += elapsed.count();
                    this->latestTimes.minITime += elapsed.count();
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
        this->printMinI();
        this->printMaxO();
        */ 
    }

    void IncrementalDFSTREEBatch::partialDFS(GraphContext *g,int w, int z) {
        if(this->marked[z] == 0 && (this->postorder[z] > this->postorder[this->u])) {
            this->marked[z] = 1;
            this->reachable[z] = 1;
            this->globalReachable[z] = 1;
            #ifndef BENCHMARK
            this->vertexCount.currentcountExploredVertices ++;
            this->vertexCount.currentTotalExploredEdges += g->getOutDegree(z);
            #endif
            for (int u : g->displayOutgoingNeighbours(z)) {
                this->partialDFS(g,z,u);
            }
            this->m = this->m+1;
            this->postorderChanged[z] = m;
     
        }
    }

    void IncrementalDFSTREEBatch::update(GraphContext *g) {
        int shift = this->m;
        int index;

        #ifndef BENCHMARK
        auto start = std::chrono::high_resolution_clock::now();
        auto end = std::chrono::high_resolution_clock::now();   
        std::chrono::duration<double> elapsed = end - start;

        start = std::chrono::high_resolution_clock::now();
        #endif
        for(int i = this->postorder[this->u]; i <= this->postorder[this->v]; i++) {
            #ifndef BENCHMARK
            this->vertexCount.currentTotalAffectedEdges += g->getOutDegree(this->reversePostorder[i]);
            this->vertexCount.currentcountAffectedVertices ++;
            #endif
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
        this->vertexCount.totalAffectedEdges += this->vertexCount.currentTotalAffectedEdges;
        #endif

        #ifndef BENCHMARK
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        this->latestTimes.currenttopologicalOrderTime += elapsed.count();
        this->latestTimes.topologicalOrderTime += elapsed.count();


        

        start = std::chrono::high_resolution_clock::now();       
        #endif
        for(int i=0; i<this->n; i++) {
            if(i == 0 || i == 1) {
                    continue;
            }
            if(this->intermediatepostorder[i]!=0) {
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
        this->vertexCount.countExploredVertices += this->vertexCount.currentcountExploredVertices;
        #endif
    }
    
    void IncrementalDFSTREEBatch::deleteReachableFromBST(GraphContext *g)
    {
        #ifndef BENCHMARK
        auto start = std::chrono::high_resolution_clock::now();
        auto end = std::chrono::high_resolution_clock::now();   
        std::chrono::duration<double> elapsed = end - start;
        #endif
        for(int i=0; i<this->n; i++) {
            if(this->globalReachable[i]!=0) {
                if(i == 0 || i == 1) {
                    continue;
                }
                #ifndef BENCHMARK
                start = std::chrono::high_resolution_clock::now();
                #endif
                for (int u : g->displayOutgoingNeighbours(i)) {                    
                       this->minIbst[u].deleteNode(i);
                }
                #ifndef BENCHMARK
                end = std::chrono::high_resolution_clock::now();
                elapsed = end - start;
                this->latestTimes.currentminITime += elapsed.count();
                this->latestTimes.minITime += elapsed.count();

                start = std::chrono::high_resolution_clock::now();
                #endif
                for (int u : g->displayIncomingNeighbours(i)) {                   
                       this->maxObst[u].deleteNode(i);
                }
                #ifndef BENCHMARK
                end = std::chrono::high_resolution_clock::now();
                elapsed = end - start;
                this->latestTimes.currentmaxOTime += elapsed.count();
                this->latestTimes.maxOTime += elapsed.count();
                #endif
            }    
        } 
    }

    void IncrementalDFSTREEBatch::insertReachableIntoBST(GraphContext *g)
    {   
        #ifndef BENCHMARK
        auto start = std::chrono::high_resolution_clock::now();
        auto end = std::chrono::high_resolution_clock::now();   
        std::chrono::duration<double> elapsed = end - start; 
        start = std::chrono::high_resolution_clock::now();
        #endif
        this->minIbst[v].addNode(u,0);
        this->minI[v] = this->minIbst[v].findMin();
        #ifndef BENCHMARK
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        this->latestTimes.currentminITime += elapsed.count();
        this->latestTimes.minITime += elapsed.count();
        
        start = std::chrono::high_resolution_clock::now();
        #endif
        this->maxObst[u].addNode(v,0);
        this->maxO[u] = this->maxObst[u].findMax(); 
        #ifndef BENCHMARK
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        this->latestTimes.currentmaxOTime += elapsed.count();
        this->latestTimes.maxOTime += elapsed.count();       
        #endif
        for(int i=0; i<this->n; i++) {
            if(this->globalReachable[i]!=0) {  
                if(i == 0 || i == 1) {
                    continue;
                } 
                #ifndef BENCHMARK
                start = std::chrono::high_resolution_clock::now();     
                #endif
                for (int u : g->displayOutgoingNeighbours(i)) {
                    if(u==1) {
                        continue;
                    }
                    this->minIbst[u].addNode(i,0);
                    this->minI[u] = (this->minIbst[u].findMin());
                }
                #ifndef BENCHMARK
                end = std::chrono::high_resolution_clock::now();
                elapsed = end - start;
                this->latestTimes.currentminITime += elapsed.count();
                this->latestTimes.minITime += elapsed.count();

                
                start = std::chrono::high_resolution_clock::now();
                #endif
                for (int u : g->displayIncomingNeighbours(i)) {
                    if (u==0) {
                        continue;
                    }
                    this->maxObst[u].addNode(i,0);
                    this->maxO[u]=  (this->maxObst[u].findMax());     
                }
                #ifndef BENCHMARK
                end = std::chrono::high_resolution_clock::now();
                elapsed = end - start;
                this->latestTimes.currentmaxOTime += elapsed.count();
                this->latestTimes.maxOTime += elapsed.count();
                #endif
            }
        }
    }

    void IncrementalDFSTREEBatch::reset() {
        std::fill(this->marked.begin(), this->marked.end(), 0);
        std::fill(this->postorderChanged.begin(), this->postorderChanged.end(), 0);
        std::fill(this->intermediatepostorder.begin(), this->intermediatepostorder.end(), 0);
        std::fill(this->reachable.begin(), this->reachable.end(), 0);
        this->m=0;
    }

    void IncrementalDFSTREEBatch::calculateMaxO(GraphContext *g) {

    }

    void IncrementalDFSTREEBatch::calculateMinI(GraphContext *g) {

    }
    void IncrementalDFSTREEBatch::printMinI() const {
        std::cout << "MinI of the graph is" << std::endl;
        for(int i=0; i < this->n; i++) {
            std::cout << i << " -> " << this->minI[i] << std::endl; 
        }
    }

    void IncrementalDFSTREEBatch::printMaxO() const {
        std::cout << "MaxO of the graph is " << std::endl;
        for(int i=0; i< this->n; i++) {
            std::cout << i << " -> " << this->maxO[i] << std::endl;
        }
    }

    void IncrementalDFSTREEBatch::printtreeminIncoming() {
        for(int i =0;i<this->n;i++) {
           std::cout << i << "->" ;
           this->minIbst[i].printTree();
           std::cout << std::endl;
        }
    }  

    void IncrementalDFSTREEBatch::printtreemaxOutgoing() {
        for(int i =0;i<this->n;i++) {
            std::cout << i << "->";
            maxObst[i].printTree();
            std::cout << std::endl;
        }
    }


    void IncrementalDFSTREEBatch::switchBSTTopo(std::vector<int>& RTopo) {
        for (int i = 0; i < this->n; i++) {
            this->minIbst[i].setTopo(RTopo);
            this->maxObst[i].setTopo(RTopo);        
        }
    }
