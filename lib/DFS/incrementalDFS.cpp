#include"incrementalDFS.h"
    
    /*
      insert function used for insertion of an edge and for maintaining all relevant datastructures like 
      reversetopological order, minI, maxO etc. It handles five cases. 
      1) edge (u,v) inserted, one new node inserted and u is new node.
      2) edge (u,v) inserted, one new node inserted and v is new node.
      3) edge (u,v) inserted, two new nodes inserted.
      4) edge (u,v) inserted, edge is from right to left.
      5) edge (u,v) inserted, edge if from left to right.  
    */
    void IncrementalDFS::insert(GraphContext *g, int u, int v) {
    	int x;
        x=this->n;
        g->insertEdge(u,v);
        this->n=g->displayNumberVertices();
        x=this->n-x;
        //std::cout << "original U= "<< u << " " << "original v= " << v << std::endl;
        this->u=g->strategy->reverseMapping[u];
        this->v=g->strategy->reverseMapping[v];
        u=g->strategy->reverseMapping[u];
        v=g->strategy->reverseMapping[v];
        //std::cout << "new U= "<< u << " " << "new v= " << v << std::endl;
        if(x>0) {
            //std::cout << " x= " << x << std::endl;
            //std::cout << " n= " << n << std::endl;
            this->reversePostorder.resize(this->n);
            this->postorder.resize(this->n);
            this->visit.resize(this->n);
            this->topOrd.resize(this->n);
            this->reverseTopologicalOrder.resize(this->n);
            this->minI.resize(this->n);
            this->maxO.resize(this->n);
            //std::cout << "after maxO ds" << std::endl;
            this->minIbst.resize(this->n);
            //std::cout << "after minIbst" << std::endl;
            this->maxObst.resize(this->n);
            //std::cout << "after maxobst" << std::endl;
            this->marked.resize(this->n);
            this->intermediatepostorder.resize(this->n);
            this->postorderChanged.resize(this->n);
            //std::cout << "after resize" << std::endl;
            if(x == 1) { 
               //std::cout << "x is equal to 1" << std::endl;
               //case where u is new node
               if(g->strategy->newu == true){
                   this->cases = 1;
                   //std::cout << "case 1" << std::endl;
                   BalancedBinaryTree binaryTreeIncoming;
                   BalancedBinaryTree binaryTreeOutgoing;
                   this->postorder[u] = this->postorder[0];
                   this->reversePostorder[postorder[u]] = u;
                   this->reverseTopologicalOrder[u] = (this->n - 1) -  this->postorder[u];
                   this->topOrd[reverseTopologicalOrder[u]] = u;           
                   binaryTreeOutgoing.addNode(this->reverseTopologicalOrder[v],v);
                   this->postorder[0] = this->postorder[u] + 1;
                   this->reversePostorder[postorder[0]] = 0;
                   this->reverseTopologicalOrder[0] = (this->n - 1) -  this->postorder[0];
                   this->topOrd[reverseTopologicalOrder[0]] = 0;
                   this->minIbst[u]= binaryTreeIncoming;
                   this->maxObst[u]= binaryTreeOutgoing;
                   this->reverseTopologicalOrder[0] = (this->n - 1) - postorder[0];
                   this->topOrd[reverseTopologicalOrder[0]] = 0;
                   this->reverseTopologicalOrder[1] = (this->n - 1) - postorder[1];
                   this->topOrd[reverseTopologicalOrder[1]] = 1;
                   
                   for(int i=2; i<this->n;i++) {
                        for (int u : g->displayOutgoingNeighbours(i)) { 
                             if(u==1) {
                                continue;
                             }                   
                             this->maxObst[i].deleteNode(reverseTopologicalOrder[u]);
                        }

                        for (int u : g->displayIncomingNeighbours(i)) {
                              if(u==0) {
                                continue;
                              }                   
                              this->minIbst[i].deleteNode(reverseTopologicalOrder[u]);
                        }
                    }    

                    for(int i=0; i<this->n; i++) {
                        this->reverseTopologicalOrder[i] = (this->n - 1) - postorder[i];
                        this->topOrd[reverseTopologicalOrder[i]] = i;
                    }
                    for(int i=2; i<this->n;i++){   
                        for (int u : g->displayOutgoingNeighbours(i)) {
                              if(u==1) {
                                continue;
                              }
                              this->maxObst[i].addNode(this->reverseTopologicalOrder[u],u);
                        }
                        for (int u : g->displayIncomingNeighbours(i)) {
                            if(u==0) {
                                continue;
                            }
                            this->minIbst[i].addNode(this->reverseTopologicalOrder[u],u);
                        }
                        if(g->strategy->incoming[i] == 2) {
                            this->minI[i] = 0;   
                        }
                        else {
                            this->minI[i] = (this->minIbst[i].findMin()).index;
                        }
                        if(g->strategy->outgoing[i] == 2) {
                            this->maxO[i] = 1;    
                        }
                        else {
                            this->maxO[i] = (this->maxObst[i].findMax()).index;
                         
                        }    
                    }
                    g->strategy->newu = false;
                }    
                else {
                    this->cases=2;
                    //case where v is new node.
                    //std::cout << "case 2" << std::endl;
                    BalancedBinaryTree binaryTreeIncoming;
                    BalancedBinaryTree binaryTreeOutgoing;
                    this->postorder[v] = this->postorder[u];
                    for(int i=2; i < this->n -1; i++) {
                       for (int u : g->displayOutgoingNeighbours(i)) {   
                          if(u == 1 || u == v) {
                            continue;
                          }                 
                          this->maxObst[i].deleteNode(this->reverseTopologicalOrder[u]);
                       }

                       for (int u : g->displayIncomingNeighbours(i)) { 
                          if(u==0) {
                            continue;
                          }                  
                          this->minIbst[i].deleteNode(this->reverseTopologicalOrder[u]);
                       }
                    }
                    for(int i=postorder[u]; i < this->n - 1; i++) {
                      this->postorder[reversePostorder[i]]++;
                    }
                    for(int i=0;i < this->n;i++) {
                      this->reversePostorder[this->postorder[i]]=i;
                      this->reverseTopologicalOrder[i]=(this->n) - 1 - this->postorder[i];
                      this->topOrd[this->reverseTopologicalOrder[i]] =i; 
                    }
                    this->minIbst[v] = binaryTreeIncoming;
                    this->maxObst[v] = binaryTreeOutgoing;
                    for(int i=2; i < this->n; i++) {
                      for (int u : g->displayOutgoingNeighbours(i)) { 
                         if(u == 1) {
                            continue;
                         }                   
                         this->maxObst[i].addNode(this->reverseTopologicalOrder[u],u);
                      }

                      for (int u : g->displayIncomingNeighbours(i)) {
                         if(u == 0) {
                            continue;
                         }                   
                         this->minIbst[i].addNode(this->reverseTopologicalOrder[u],u);
                      }
                      
                      if(g->strategy->incoming[i] == 2) {
                            this->minI[i] = 0;   
                      }
                      else {
                            this->minI[i] = (this->minIbst[i].findMin()).index;
                      }
                      if(g->strategy->outgoing[i] == 2) {
                            this->maxO[i] = 1;    
                      }
                      else {
                            this->maxO[i] = (this->maxObst[i].findMax()).index;
                         
                      }   
                    }
                    g->strategy->newv= false;    
                }    
            }
            else {
               this->cases = 3;
               //std::cout << "if x is 2" << std::endl; 
               BalancedBinaryTree binaryTreeIncomingu;
               BalancedBinaryTree binaryTreeOutgoingu;
               BalancedBinaryTree binaryTreeIncomingv;
               BalancedBinaryTree binaryTreeOutgoingv;


               this->postorder[v] = this->postorder[0];
               this->reversePostorder[postorder[v]] = v;
               this->reverseTopologicalOrder[v] = (this->n - 1) -  this->postorder[v];
               this->topOrd[reverseTopologicalOrder[v]] = v;
               binaryTreeOutgoingu.addNode(this->reverseTopologicalOrder[v],v);
         
               this->postorder[u] = this->postorder[v] + 1;
               this->reversePostorder[postorder[u]] = u;
               this->reverseTopologicalOrder[u] = (this->n - 1) -  this->postorder[u];
               this->topOrd[reverseTopologicalOrder[u]] = u;          
               binaryTreeIncomingv.addNode(this->reverseTopologicalOrder[u],u);

               this->postorder[0] = this->postorder[u] + 1;
               this->reversePostorder[postorder[0]] = 0;
               this->reverseTopologicalOrder[0] = (this->n - 1) -  this->postorder[0];
               this->topOrd[reverseTopologicalOrder[0]] = 0;
               binaryTreeOutgoingu.addNode(this->reverseTopologicalOrder[v],v);
               binaryTreeIncomingv.addNode(this->reverseTopologicalOrder[u],u);
               this->minIbst[u]= binaryTreeIncomingu;
               this->minIbst[v]= binaryTreeIncomingv;
               this->maxObst[u]= binaryTreeOutgoingu;
               this->maxObst[v]= binaryTreeOutgoingv;
               for(int i=2; i<this->n - 2;i++) {
                   for (int u : g->displayOutgoingNeighbours(i)) {  
                       if(u == 1) {
                         continue;
                       }                  
                       this->maxObst[i].deleteNode(this->reverseTopologicalOrder[u]);
                   }

                   for (int u : g->displayIncomingNeighbours(i)) {  
                       if(u == 0) {
                         continue;
                       }                 
                       this->minIbst[i].deleteNode(this->reverseTopologicalOrder[u]);
                   }
                }   
                for(int i=0; i<this->n;i++) {

                   this->reverseTopologicalOrder[i] = (this->n - 1) - postorder[i];
                   this->topOrd[reverseTopologicalOrder[i]] = i;
                }

                for(int i=2; i<this->n; i++) {    
                   for (int u : g->displayOutgoingNeighbours(i)) {
                      if(u==1) {
                        continue;
                      }
                      this->maxObst[i].addNode(this->reverseTopologicalOrder[u],u);
                   }
                   
                   for (int u : g->displayIncomingNeighbours(i)) {
                      if(u==0) {
                        continue;
                      }
                      this->minIbst[i].addNode(this->reverseTopologicalOrder[u],u);
                   } 
                   if(g->strategy->incoming[i] == 2) {
                        this->minI[i] = 0;   
                   }
                   else {
                        this->minI[i] = (this->minIbst[i].findMin()).index;
                   }
                   if(g->strategy->outgoing[i] == 2) {
                        this->maxO[i] = 1;    
                   }
                   else {
                        this->maxO[i] = (this->maxObst[i].findMax()).index;
                         
                   }   
                }
                g->strategy->newu = false;
                g->strategy->newv = false;
                           
            }
        }    
        else {
            this->cases = 4;
            //std::cout << "if no new node added" << std::endl;
            if(this->postorder[u] < this->postorder[v]) {
    		    this->partialDFS(g,u,v);
    		    this->update(g);
                this->reset();
    	    }
            else {
                this->maxObst[u].addNode(this->reverseTopologicalOrder[v],v);
                this->minIbst[v].addNode(this->reverseTopologicalOrder[u],u);
                this->maxO[u] = (this->maxObst[u].findMax()).index;
                this->minI[v] = (this->minIbst[v].findMin()).index;
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

    void IncrementalDFS::partialDFS(GraphContext *g,int w, int z) {
    	if(this->marked[z] == 0 && (this->postorder[z] > this->postorder[this->u])) {
            this->marked[z] =1;
            for (int u : g->displayOutgoingNeighbours(z)) {
                this->partialDFS(g,z,u);
            }
            this->m = this->m+1;
            this->postorderChanged[z] = m;
     
    	}
    }

    void IncrementalDFS::update(GraphContext *g) {
        int shift = this->m;
    	int index;
    	for(int i = this->postorder[this->u]; i <= this->postorder[this->v]; i++) {
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
    	for(int i=0; i<this->n; i++) {
            if(this->intermediatepostorder[i]!=0) {
                //std::cout << "i = " << i << std::endl;
                if(i == 0 || i == 1) {
                    continue;
                }
                for (int u : g->displayOutgoingNeighbours(i)) {
                       if(u == 1) {
                         continue;
                       }                    
                       this->minIbst[u].deleteNode(this->reverseTopologicalOrder[i]);
                }

                for (int u : g->displayIncomingNeighbours(i)) {
                       if( u == 0) {
                          continue;
                       }                   
                       this->maxObst[u].deleteNode(this->reverseTopologicalOrder[i]);
                }
            }    
        } 
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
        for(int i=0; i<this->n; i++) {
            if(this->intermediatepostorder[i]!=0) {
                if(i == 0 || i == 1) {
                    continue;
                }        
                for (int u : g->displayOutgoingNeighbours(i)) {
                    if(u == 1) {
                        continue;
                    }
                    this->minIbst[u].addNode(this->reverseTopologicalOrder[i],i);
                    this->minI[u] = (this->minIbst[u].findMin()).index;
                    
                }

                for (int u : g->displayIncomingNeighbours(i)) {
                    if( u == 0) {
                        continue;
                    }
                    this->maxObst[u].addNode(this->reverseTopologicalOrder[i],i);  
                    this->maxO[u] = (this->maxObst[u].findMax()).index;
                }
            }
    	}

    }
    
    void IncrementalDFS::reset() {
        std::fill(this->marked.begin(), this->marked.end(), 0);
        std::fill(this->postorderChanged.begin(), this->postorderChanged.end(), 0);
        std::fill(this->intermediatepostorder.begin(), this->intermediatepostorder.end(), 0);
        this->m=0;
    }

    void IncrementalDFS::calculateMaxO(GraphContext *g) {

    }

    void IncrementalDFS::calculateMinI(GraphContext *g) {

    }
    void IncrementalDFS::printMinI() const {
        std::cout << "MinI of the graph is" << std::endl;
        for(int i=0; i < this->n; i++) {
            std::cout << i << " -> " << this->minI[i] << std::endl; 
        }
    }

    void IncrementalDFS::printMaxO() const {
        std::cout << "MaxO of the graph is " << std::endl;
        for(int i=0; i< this->n; i++) {
            std::cout << i << " -> " << this->maxO[i] << std::endl;
        }
    }

    void IncrementalDFS::printtreeminIncoming() {
        for(int i =0;i<this->n;i++) {
           std::cout << i << "->" ;
           this->minIbst[i].printTree();
           std::cout << std::endl;
        }
    }  

    void IncrementalDFS::printtreemaxOutgoing() {
        for(int i =0;i<this->n;i++) {
            std::cout << i << "->";
            maxObst[i].printTree();
            std::cout << std::endl;
        }
    }