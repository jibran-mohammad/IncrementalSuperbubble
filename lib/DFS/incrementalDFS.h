#ifndef INCREMENTAL_DFS
#define INCREMENTAL_DFS
#include "baseDFS.h"
#include "minImaxOdatastructure.h"
#include"../tree/differentAVLTree.h"
#include <climits>
#include <chrono>
/*
  Incremental DFS algorithm, it uses incremental DFS in DAG algorithm and uses different balanced binary search
  tree for storing incoming and outgoing edges according to the topological order. Theoritically it takes 
  O(m^2logn) for insertion of 'm' edges but practically it takes less time.
*/
class IncrementalDFS:public BaseDFS , public IDatastructure  {
public:
	  struct IncrementalUpdateTimes {
    double insertTime = 0.0;
    double dfsTime = 0.0;
    double topologicalOrderTime = 0.0;
    double minITime = 0.0;
    double maxOTime = 0.0;
    double currentinsertTime = 0.0;
    double currentdfsTime = 0.0;
    double currenttopologicalOrderTime = 0.0;
    double currentminITime = 0.0;
    double currentmaxOTime = 0.0;

    };
    struct TypeOfEdgeInsertion {
    	int case1 = 0;
    	int case2 = 0;
    	int case3 = 0;
    	int case4 = 0;
    	int case5 = 0;
    	int currentcase1 = 0;
    	int currentcase2 = 0;
    	int currentcase3 = 0;
    	int currentcase4 = 0;
    	int currentcase5 = 0;
    };
    IncrementalUpdateTimes latestTimes;
    TypeOfEdgeInsertion edgeTypes;
    std::vector<int> marked;
    std::vector<int> postorderChanged;
    std::vector<int> intermediatepostorder;
    std::vector<DifferentBalancedBinaryTree> minIbst;
	  std::vector<DifferentBalancedBinaryTree> maxObst;
    int m,u,v;    

public:
	IncrementalDFS(GraphContext *g): BaseDFS(g) {
        cases = 0;
        g->insertEdge(0,1);
        n=g->displayNumberVertices();
        this->m=0;
        DifferentBalancedBinaryTree binaryTreeIncoming1;
        DifferentBalancedBinaryTree binaryTreeOutgoing1;
        DifferentBalancedBinaryTree binaryTreeIncoming2;
        DifferentBalancedBinaryTree binaryTreeOutgoing2;
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