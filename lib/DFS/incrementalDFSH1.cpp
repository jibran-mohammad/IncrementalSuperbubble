#include"incrementalDFSAVLTree.h"
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
    void IncrementalDFSTREE::insert(GraphContext *g, int u, int v) {
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
            this->reachable.resize(this->n);
            this->minI.resize(this->n);
            this->maxO.resize(this->n);
            this->minIbst.resize(this->n);
            this->maxObst.resize(this->n);
            this->marked.resize(this->n);
            this->intermediatepostorder.resize(this->n);
            this->postorderChanged.resize(this->n);
            if(x == 1) {
               //std::cout << "x is equal to 1" << std::endl;
               //case where u is new node
                this->cases = 1;
               if(g->strategy->newu == true){
                   //std::cout << "case 1" << std::endl;
                   DifferentBalancedBinaryTree binaryTreeIncoming(this->reverseTopologicalOrder);
                   DifferentBalancedBinaryTree binaryTreeOutgoing(this->reverseTopologicalOrder);
                   this->postorder[u] = this->postorder[0];
                   this->reversePostorder[postorder[u]] = u;
                   this->reverseTopologicalOrder[u] = (this->n - 1) -  this->postorder[u];
                   this->topOrd[reverseTopologicalOrder[u]] = u;           
                   
                   this->postorder[0] = this->postorder[u] + 1;
                   this->reversePostorder[postorder[0]] = 0;
                   this->reverseTopologicalOrder[0] = (this->n - 1) -  this->postorder[0];
                   this->topOrd[reverseTopologicalOrder[0]] = 0;

                   this->minIbst[u]= binaryTreeIncoming;
                   this->maxObst[u]= binaryTreeOutgoing;
                
                   for(int i=0; i<this->n;i++) {
                      this->reverseTopologicalOrder[i] = (this->n - 1) - postorder[i];
                      this->topOrd[reverseTopologicalOrder[i]] = i;
                   }
                   this->maxObst[u].addNode(v,0);
                   this->maxO[u] =  this->maxObst[u].findMax();
                   this->minIbst[v].addNode(u,0);
                   this->minI[v] = this->minIbst[v].findMin();
                   g->strategy->newu = false;
               }    
               else {
                    //v is the new node
                    //std::cout << "case 2" << std::endl;
                    this->cases=2;
                    DifferentBalancedBinaryTree binaryTreeIncoming(this->reverseTopologicalOrder);
                    DifferentBalancedBinaryTree binaryTreeOutgoing(this->reverseTopologicalOrder);
                    this->postorder[v] = this->postorder[u];
                    for(int i=postorder[u]; i < this->n - 1; i++) {
                      this->postorder[reversePostorder[i]]++;
                    }
                    for(int i=0;i < this->n;i++) {
                      this->reversePostorder[this->postorder[i]]=i;
                      this->reverseTopologicalOrder[i]=(this->n) - 1 - this->postorder[i];
                      this->topOrd[this->reverseTopologicalOrder[i]] =i; 
                    }
                    binaryTreeIncoming.addNode(u,0);
                    this->maxObst[u].addNode(v,0);
                    this->minIbst[v] = binaryTreeIncoming;
                    this->maxObst[v] = binaryTreeOutgoing;
                    this->maxO[u] = this->maxObst[u].findMax();
                    this->minI[v] = this->minIbst[v].findMin();
                    this->maxO[v] = 1;
                    g->strategy->newv= false;
                }
            }
            else {
               //std::cout << "if x is 2" << std::endl; 
               this->cases = 3;
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
               DifferentBalancedBinaryTree binaryTreeIncomingu(this->reverseTopologicalOrder);
               DifferentBalancedBinaryTree binaryTreeOutgoingu(this->reverseTopologicalOrder);
               DifferentBalancedBinaryTree binaryTreeIncomingv(this->reverseTopologicalOrder);
               DifferentBalancedBinaryTree binaryTreeOutgoingv(this->reverseTopologicalOrder);
               binaryTreeOutgoingu.addNode(v,0);
               binaryTreeIncomingv.addNode(u,0);
               this->minIbst[u]= binaryTreeIncomingu;
               this->minIbst[v]= binaryTreeIncomingv;
               this->maxObst[u]= binaryTreeOutgoingu;
               this->maxObst[v]= binaryTreeOutgoingv;
               this->minI[u] = 0;
               this->maxO[u] = v;
               this->minI[v] = u;
               this->maxO[v] = 1; 
               g->strategy->newu = false;
               g->strategy->newv = false;
               }           
            }
        else {
            this->cases = 4;
            //std::cout << "if no new node added" << std::endl;
            if(this->postorder[u] < this->postorder[v]) {
                //std::cout << "left to right edge" << std::endl;
                this->partialDFS(g,u,v);
                this->update(g);
                this->reset();
            }
            else {
                //std::cout << "right to left  edge" << std::endl;
                this->maxObst[u].addNode(v,0);
                this->minIbst[v].addNode(u,0);
                this->maxO[u] = (this->maxObst[u].findMax());
                this->minI[v] = (this->minIbst[v].findMin());
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

    void IncrementalDFSTREE::partialDFS(GraphContext *g,int w, int z) {
        if(this->marked[z] == 0 && (this->postorder[z] > this->postorder[this->u])) {
            this->marked[z] = 1;
            this->reachable[z] = 1;
            for (int u : g->displayOutgoingNeighbours(z)) {
                this->partialDFS(g,z,u);
            }
            this->m = this->m+1;
            this->postorderChanged[z] = m;
     
        }
    }

    void IncrementalDFSTREE::update(GraphContext *g) {
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
            if(this->reachable[i]!=0) {
                if(i == 0 || i == 1) {
                    continue;
                }
                for (int u : g->displayOutgoingNeighbours(i)) {                    
                       this->minIbst[u].deleteNode(i);
                }
            }    
        }        
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
        this->minIbst[v].addNode(u,0);
        this->maxObst[u].addNode(v,0);
        this->minI[v] = this->minIbst[v].findMin();
        this->maxO[u] = this->maxObst[u].findMax();        
        for(int i=0; i<this->n; i++) {
            if(this->reachable[i]!=0) {  
                if(i == 0 || i == 1) {
                    continue;
                }      
                for (int u : g->displayOutgoingNeighbours(i)) {
                    if(u==1) {
                        continue;
                    }
                    this->minIbst[u].addNode(i,0);
                    this->minI[u] = (this->minIbst[u].findMin());    
                }
                for(int u: g->displayIncomingNeighbours(i)) {
                    if(u==0) {
                        continue;
                    }
                    if(this->reverseTopologicalOrder[this->maxO[u]] < this->reverseTopologicalOrder[i] ) {
                        this->maxO[u] = i;
                    }
                }
            }
        }

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
