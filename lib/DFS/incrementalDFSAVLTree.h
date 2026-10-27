#ifndef INCREMENTAL_DFS_AVL_TREE
#define INCREMENTAL_DFS_AVL_TREE
#include "baseDFS.h"
#include "minImaxOdatastructure.h"
#include"../tree/differentAVLTree.h"
#include <climits>
#include <chrono>
/*
 class IncremetalDFSTREE, it uses incremental DFS in DAG algorithm, and a special type of balanced binary 
 search tree which stores the indices of the graph as key, and has access to reverse topological order.This
 algorithm takes O(mnlogn) time.
*/
class IncrementalDFSTREE:public BaseDFS , public IDatastructure  {
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
    struct CountOfVertices {   // 64 bit: perBatch sums over a whole HPRC genome exceed 2^31
    	long long countAffectedVertices = 0;
    	long long countExploredVertices = 0;
    	long long currentcountAffectedVertices = 0;
    	long long currentcountExploredVertices = 0;
    	long long totalEdges = 0;
    	long long currentTotalEdges =0;
    	long long totalExploredEdges = 0;
    	long long totalAffectedEdges = 0;
    	long long currentTotalExploredEdges =0;
    	long long currentTotalAffectedEdges =0;
    };
    IncrementalUpdateTimes latestTimes;
    TypeOfEdgeInsertion edgeTypes;
    CountOfVertices vertexCount;
    std::vector<int> marked;
    std::vector<int> postorderChanged;
    std::vector<int> intermediatepostorder;
    std::vector<int> reachable;
    std::vector<DifferentBalancedBinaryTree> minIbst;
	std::vector<DifferentBalancedBinaryTree> maxObst;
    int m,u,v;    
    // number of updates after which the whole order was recomputed with the
    // static DFS (rebuildFromStaticDFS)
    int staticRebuilds;

public:
	IncrementalDFSTREE(GraphContext *g): BaseDFS(g) {
        cases=0;
        g->insertEdge(0,1);
        n=g->displayNumberVertices();
        this->m=0;
        this->staticRebuilds=0;
    	this->reverseTopologicalOrder.resize(n);
    	this->topOrd.resize(n);
    	this->minI.resize(n);
		this->maxO.resize(n);
		this->reachable.resize(n);
		this->minIbst.resize(n);
		this->maxObst.resize(n);
		this->reverseTopologicalOrder[0] = 0;
		this->topOrd[0] = 0;
		this->reverseTopologicalOrder[1] = 1;
		this->topOrd[1] = 1;
		this->minI[0] = INT_MAX;
		this->minI[1] = 0;
		this->maxO[0] = 1;
		this->maxO[1] = INT_MAX;
		DifferentBalancedBinaryTree binaryTreeIncoming1(this->reverseTopologicalOrder);
        DifferentBalancedBinaryTree binaryTreeOutgoing1(this->reverseTopologicalOrder);
        DifferentBalancedBinaryTree binaryTreeIncoming2(this->reverseTopologicalOrder);
        DifferentBalancedBinaryTree binaryTreeOutgoing2(this->reverseTopologicalOrder);
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
	void rebuildFromStaticDFS(GraphContext *g);
	void calculateMinI(GraphContext *g) override;
	void calculateMaxO(GraphContext *g) override;
	void printMinI() const override;
	void printMaxO() const override;
	void printtreeminIncoming();
	void printtreemaxOutgoing();
	void deleteReachableFromBST(GraphContext *g);
	void insertReachableIntoBST(GraphContext *g);

};
#endif