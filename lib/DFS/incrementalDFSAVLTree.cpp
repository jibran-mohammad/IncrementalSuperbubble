#include"incrementalDFSAVLTree.h"
    /*
      insert function used for insertion of an edge and for maintaining all relevant datastructures like 
      reversetopological order, minI, maxO etc. It handles five cases. 
      1) edge (u,v) inserted, one new node inserted and u is new node.
      2) edge (u,v) inserted, one new node inserted and v is new node.
      3) edge (u,v) inserted, two new nodes inserted.
      4) edge (u,v) inserted, both nodes exist, edge is from left to right in the DFS tree
         (u is visited before v's part of the tree, u after v in the topological order:
         the order is broken and partialDFS re-places v and what it reaches after u).
      5) edge (u,v) inserted, both nodes exist, edge is from right to left in the DFS tree
         (u before v in the topological order: the order stays valid, nothing is reordered).

      Case 1: u is placed as the rightmost child of 0 (directly after 0 in
      the topological order), so (u,v) is a right-to-left edge and nothing is
      reordered: all other positions move up by one (O(n)), no partialDFS,
      the minI trees only get u in v's tree and u gets a maxO tree with v. If
      v's only parent was 0, u takes v's place instead (directly before v),
      again without a rebuild. (The same as in incMinIBST.)

      Case 5 when v's only parent was the artificial source 0: the graph
      deletes the edge 0->v (g->strategy->artificialParentRemoved), but in the
      DFS tree v and its subtree are still children of 0, so they are not placed
      under u: superbubbles that now contain v are not contiguous in the order
      and computeSB misses them. The order is then recomputed once with the
      static DFS (rebuildFromStaticDFS), in which v's subtree hangs under u.
      the difference from the O(m^2logn) is the number of vertices we need to explore whose topological order
      got changed due to insertion of an edge. We are only changing the vertices in the special balanced binary
      search tree whos reachability changed.  
    */
    void IncrementalDFSTREE::insert(GraphContext *g, int u, int v) {
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
        // perBatch: the graph size after the genome (like TotalVertices = n),
        // not summed over the genome's edges (affected / explored are summed:
        // that is the work of the genome)
        this->vertexCount.totalEdges = this->vertexCount.currentTotalEdges;
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
            this->visit.resize(this->n);
            this->topOrd.resize(this->n);
            this->reverseTopologicalOrder.resize(this->n);
            this->reachable.resize(this->n);
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
                   {
                   // v has another parent: u becomes the rightmost child of 0
                   // in the DFS tree (visited last), so it comes directly after
                   // 0 in the topological order and (u,v) is a right-to-left
                   // edge: nothing else is reordered, no partialDFS.
                   // v's only parent was 0 (the edge 0->v is deleted): u takes
                   // v's place as a child of 0, i.e. directly before v; v's
                   // subtree is visited from u at the same point as before
                   // from 0, so it stays the same and no rebuild is needed.
                   // Every vertex from that position on moves one position
                   // later; their relative order is unchanged, so the minI and
                   // maxO trees (ordered by the current positions) stay valid.
                   // (As in incMinIBST, lib/DFS/incrementalDFSAVLTreeH4.cpp.)
                   int pos = g->strategy->artificialParentRemoved ? this->reverseTopologicalOrder[v] : 1;
                   for(int i = this->n - 2; i >= pos; i--) {
                       this->topOrd[i + 1] = this->topOrd[i];
                       this->reverseTopologicalOrder[this->topOrd[i + 1]] = i + 1;
                   }
                   this->reverseTopologicalOrder[u] = pos;
                   this->topOrd[pos] = u;
                   #ifndef BENCHMARK
                   end = std::chrono::high_resolution_clock::now();
                   elapsed = end - start;
                   this->latestTimes.currenttopologicalOrderTime += elapsed.count();
                   this->latestTimes.topologicalOrderTime += elapsed.count();

                   start = std::chrono::high_resolution_clock::now();
                   #endif
                   // u's only child is v
                   this->maxObst[u] = DifferentBalancedBinaryTree(this->reverseTopologicalOrder);
                   this->maxObst[u].addNode(v,0);
                   this->maxO[u] = this->maxObst[u].findMax();
                   #ifndef BENCHMARK
                   end = std::chrono::high_resolution_clock::now();
                   elapsed = end - start;
                   this->latestTimes.currentmaxOTime += elapsed.count();
                   this->latestTimes.maxOTime += elapsed.count();

                   start = std::chrono::high_resolution_clock::now();
                   #endif
                   // u has no parent but 0 (minI[u] = 0); v gets the parent u,
                   // which is now its earliest parent
                   this->minI[u] = 0;
                   this->minIbst[u] = DifferentBalancedBinaryTree(this->reverseTopologicalOrder);
                   this->minIbst[v].addNode(u,0);
                   this->minI[v] = this->minIbst[v].findMin();
                   #ifndef BENCHMARK
                   end = std::chrono::high_resolution_clock::now();
                   elapsed = end - start;
                   this->latestTimes.currentminITime += elapsed.count();
                   this->latestTimes.minITime += elapsed.count();
                   #endif
                   }
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
                    int i =0;
                    for( i=this->n-2; i>this->reverseTopologicalOrder[this->u];i--) {
                        this->reverseTopologicalOrder[this->topOrd[i]]++;
                    }
                    this->reverseTopologicalOrder[v] = this->reverseTopologicalOrder[this->u] + 1;
                    i=0;
                    for(i=this->n-1; i>this->reverseTopologicalOrder[this->u] + 1;i--) {
                        this->topOrd[i] = this->topOrd[i-1];

                    }
                    this->topOrd[i] = this->v;
                    i=0;
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
               this->reverseTopologicalOrder[1] = this->n - 1;
               this->topOrd[this->reverseTopologicalOrder[1]] = 1;

               this->reverseTopologicalOrder[v] = this->reverseTopologicalOrder[1] - 1;
               this->topOrd[this->reverseTopologicalOrder[v]] = v;
            
               this->reverseTopologicalOrder[u] = this->reverseTopologicalOrder[1] - 2;
               this->topOrd[reverseTopologicalOrder[u]] = u;  
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
            if(this->reverseTopologicalOrder[u] > this->reverseTopologicalOrder[v]) {
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
                this->deleteReachableFromBST(g);

                this->update(g);

                this->insertReachableIntoBST(g);
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
                if(g->strategy->artificialParentRemoved) {
                    // v's only parent was 0: move v's subtree under u
                    this->rebuildFromStaticDFS(g);
                }
                else {
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
        }
        
    }

    void IncrementalDFSTREE::partialDFS(GraphContext *g,int w, int z) {
        //std::cout << "inside partial DFS" << std::endl;
        //std::cout << marked[z] << std::endl;
        if(this->marked[z] == 0 && (this->reverseTopologicalOrder[z] < this->reverseTopologicalOrder[this->u])) {
            //std::cout << z << std::endl;
            this->marked[z] =1;
            this->reachable[z] = 1;
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

    void IncrementalDFSTREE::update(GraphContext *g) {
        int shift = this->m;
        //std::cout << "shift = " << shift << std::endl;
        int index;
        
        #ifndef BENCHMARK
        auto start = std::chrono::high_resolution_clock::now();
        auto end = std::chrono::high_resolution_clock::now();   
        std::chrono::duration<double> elapsed = end - start;


        start = std::chrono::high_resolution_clock::now();
        #endif
        
        for(int i = this->reverseTopologicalOrder[this->u]; i >= this->reverseTopologicalOrder[this->v]; i--) {
            #ifndef BENCHMARK
            // i is a position: the vertex there is topOrd[i] (reverseTopologicalOrder
            // maps a vertex to its position; it used to be read with i, i.e. the
            // out-degree of an arbitrary vertex was added)
            this->vertexCount.currentTotalAffectedEdges += g->getOutDegree(this->topOrd[i]);
            this->vertexCount.currentcountAffectedVertices ++;
            #endif
            index = this->topOrd[i];
            
            if(this->marked[index]==1) {
                this->intermediatepostorder[index] = this->reverseTopologicalOrder[this->u] - this->postorderChanged[index] + 1;
                this->marked[index] = 0;
            }
            else {
                this->intermediatepostorder[index] = this->reverseTopologicalOrder[this->u] - shift;
                shift = shift + 1;
            }
            //std::cout << index << " " << intermediatepostorder[index] << std::endl ;
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
            if(this->intermediatepostorder[i]!=0) {
                //std::cout << "i = " << i << std::endl;
                if(i == 0) {
                    continue;
                }       
                this->reverseTopologicalOrder[i]= this->intermediatepostorder[i];
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

    
    
    void IncrementalDFSTREE::deleteReachableFromBST(GraphContext *g)
    {
        #ifndef BENCHMARK
        auto start = std::chrono::high_resolution_clock::now();
        auto end = std::chrono::high_resolution_clock::now();   
        std::chrono::duration<double> elapsed = end - start;
        #endif
        for(int i=0; i<this->n; i++) {
            if(this->reachable[i]!=0) {
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

    void IncrementalDFSTREE::insertReachableIntoBST(GraphContext *g)
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
            if(this->reachable[i]!=0) {  
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

    /*
      Recomputes the whole order with the static DFS (BaseDFS::dfs) on the
      current graph, then all minI and maxO trees and minI/maxO from it. Used in cases 1 and 5 when the
      inserted edge removed v's artificial edge 0->v (see insert).
    */
    void IncrementalDFSTREE::rebuildFromStaticDFS(GraphContext *g) {
        #ifndef BENCHMARK
        auto start = std::chrono::high_resolution_clock::now();
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        start = std::chrono::high_resolution_clock::now();
        #endif
        std::fill(this->visit.begin(), this->visit.end(), 0);
        this->postorderNumber = 0;
        this->dfs(g);
        #ifndef BENCHMARK
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        this->latestTimes.currentdfsTime += elapsed.count();
        this->latestTimes.dfsTime += elapsed.count();
        start = std::chrono::high_resolution_clock::now();
        #endif
        // the trees are ordered by the old order: rebuilt from the real edges
        // (the artificial vertices 0 and 1 are never stored in them)
        for(int i = 2; i < this->n; i++) {
            DifferentBalancedBinaryTree incoming(this->reverseTopologicalOrder);
            for(int p : g->displayIncomingNeighbours(i)) {
                if(p != 0) {
                    incoming.addNode(p,0);
                }
            }
            this->minIbst[i] = incoming;
            this->minI[i] = this->minIbst[i].findMin();     // 0 if no real parent
        }
        #ifndef BENCHMARK
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        this->latestTimes.currentminITime += elapsed.count();
        this->latestTimes.minITime += elapsed.count();
        start = std::chrono::high_resolution_clock::now();
        #endif
        for(int i = 2; i < this->n; i++) {
            DifferentBalancedBinaryTree outgoing(this->reverseTopologicalOrder);
            for(int c : g->displayOutgoingNeighbours(i)) {
                if(c != 1) {
                    outgoing.addNode(c,0);
                }
            }
            this->maxObst[i] = outgoing;
            this->maxO[i] = this->maxObst[i].findMax();     // 1 if no real child
        }
        #ifndef BENCHMARK
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        this->latestTimes.currentmaxOTime += elapsed.count();
        this->latestTimes.maxOTime += elapsed.count();
        #endif
        this->reset();
        this->staticRebuilds++;
        #ifdef DEBUG_OUTPUT  // debug output: compiled only with `make ... DEBUG=1`
        std::cout << "v's only parent was 0: order recomputed with the static DFS" << std::endl;
        #endif
    }

    void IncrementalDFSTREE::reset() {
        std::fill(this->marked.begin(), this->marked.end(), 0);
        std::fill(this->postorderChanged.begin(), this->postorderChanged.end(), 0);
        std::fill(this->intermediatepostorder.begin(), this->intermediatepostorder.end(), 0);
        std::fill(this->reachable.begin(), this->reachable.end(), 0);
        this->m=0;
    }

    void IncrementalDFSTREE::calculateMaxO(GraphContext *g) {

    }

    void IncrementalDFSTREE::calculateMinI(GraphContext *g) {

    }
    void IncrementalDFSTREE::printMinI() const {
        std::cout << "MinI of the graph is" << std::endl;
        for(int i=0; i < this->n; i++) {
            std::cout << i << " -> " << this->minI[i] << std::endl; 
        }
    }

    void IncrementalDFSTREE::printMaxO() const {
        std::cout << "MaxO of the graph is " << std::endl;
        for(int i=0; i< this->n; i++) {
            std::cout << i << " -> " << this->maxO[i] << std::endl;
        }
    }

    void IncrementalDFSTREE::printtreeminIncoming() {
        for(int i =0;i<this->n;i++) {
           std::cout << i << "->" ;
           this->minIbst[i].printTree();
           std::cout << std::endl;
        }
    }  

    void IncrementalDFSTREE::printtreemaxOutgoing() {
        for(int i =0;i<this->n;i++) {
            std::cout << i << "->";
            maxObst[i].printTree();
            std::cout << std::endl;
        }
    }
