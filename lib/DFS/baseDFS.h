#ifndef BASEDFS_
#define BASEDFS_
#include "dfs.h"
#include <chrono>


class BaseDFS: public IDFS {
public:
	//int a,b;
	BaseDFS(GraphContext *g) {
    	n=g->displayNumberVertices();
    	this->postorderNumber = 0;
    	this->postorder.resize(n);
    	this->reverseTopologicalOrder.resize(n);
    	this->visit.resize(n);
    	this->topOrd.resize(n);
    	this->reversePostorder.resize(n);
    	this->minI.resize(n);
		this->maxO.resize(n);
		this->visit.assign(this->n, 0); 
        this->reversePostorder.clear();
	}
	struct IncrementalUpdateTimes {
	    double dfsTime = 0.0;
	    double topologicalOrderTime = 0.0;
	    double currentdfsTime = 0.0;
	    double currenttopologicalOrderTime = 0.0;
    };
    IncrementalUpdateTimes latestTimes;
public:
	void dfs(GraphContext *g) override;
	void dfsUtil(GraphContext *G, int v);
	void printTopologicalOrder() const override;
	void printReverseTopologicalOrder() const override;
	void printPostOrder() const override;
	void insert(GraphContext * g, int u , int v) override;
	void printreversePostOrder() const override;
	void reset() override;
};
#endif