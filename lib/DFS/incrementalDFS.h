#ifndef INCREMENTAL_DFS
#define INCREMENTAL_DFS
#include "baseDFS.h"
#include "minImaxOdatastructure.h"
#include"../tree/avlTree.h"
#include <climits>
/*
  Incremental DFS algorithm, it uses incremental DFS in DAG algorithm and uses simple balanced binary search
  tree for storing incoming and outgoing edges according to the topological order. Theoritically it takes 
  O(m^2logn) for insertion of 'm' edges but practically it takes less time.
*/
class IncrementalDFS:public BaseDFS , public IDatastructure  {
public:
    std::vector<int> marked;
    std::vector<int> postorderChanged;
    std::vector<int> intermediatepostorder;
    std::vector<BalancedBinaryTree> minIbst;
	std::vector<BalancedBinaryTree> maxObst;
    int m,u,v;    

public:
	IncrementalDFS(GraphContext *g): BaseDFS(g) {
        cases = 0;
        g->insertEdge(0,1);
        n=g->displayNumberVertices();
        this->m=0;
        BalancedBinaryTree binaryTreeIncoming1;
        BalancedBinaryTree binaryTreeOutgoing1;
        BalancedBinaryTree binaryTreeIncoming2;
        BalancedBinaryTree binaryTreeOutgoing2;
		this->postorder.resize(n);
    	this->reverseTopologicalOrder.resize(n);
    	this->topOrd.resize(n);
    	this->reversePostorder.resize(n);
    	this->minI.resize(n);
		this->maxO.resize(n);
		this->minIbst.resize(n);
		this->maxObst.resize(n);
		this->postorder[0] = 1;
		this->reversePostorder[0] = 1;
		this->reverseTopologicalOrder[0] = 0;
		this->topOrd[0] = 0;
		this->postorder[1] = 0;
		this->reversePostorder[1] = 0;
		this->reverseTopologicalOrder[1] = 1;
		this->topOrd[1] = 1;
		this->minI[0] = INT_MAX;
		this->minI[1] = 0;
		this->maxO[0] = 1;
		this->maxO[1] = INT_MAX;
	    this->minIbst[0] = binaryTreeIncoming1;
	    this->minIbst[1] = binaryTreeIncoming2;
	    this->maxObst[0] = binaryTreeOutgoing1;
	    this->maxObst[1] = binaryTreeOutgoing2;
	    //this->printtreeminIncoming();
	    //this->printtreemaxOutgoing();
		marked.resize(n);
		postorderChanged.resize(n);
		intermediatepostorder.resize(n);
		for(int i=0; i<n; i++) {
			marked[i] = 0;
		}


    }
	void insert(GraphContext *g,int u, int v) override;
	void reset() override;
	void partialDFS(GraphContext *g,int u , int v );
	void update(GraphContext *g);
	void calculateMinI(GraphContext *g) override;
	void calculateMaxO(GraphContext *g) override;
	void printMinI() const override;
	void printMaxO() const override;
	void printtreeminIncoming();
	void printtreemaxOutgoing();

};
#endif