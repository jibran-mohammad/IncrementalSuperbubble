#ifndef CONTEXTDFS_
#define CONTEXTDFS_
#include "simpleDFS.h"
#include "../graph/contextGraph.h"
/*
 DFSContext class used for performing all the operations provided by all the different types of DFS algorithms
 someone who needs to use any algorithm has to go through this class to access, whatever functionality is 
 provided.It takes object of the algorithm as the argument and provides necessary functionality. 
*/
class DFSCONTEXT {
public: 
	IDFS *dfsStrategy;
	IDatastructure *datastructureStrategy;
public:
	DFSCONTEXT(IDFS *strategy):dfsStrategy(strategy), datastructureStrategy(dynamic_cast<IDatastructure*>(strategy)) {}

	void setStrategy(IDFS *strategy) {
		this->dfsStrategy = strategy;
		this->datastructureStrategy = dynamic_cast<IDatastructure*>(strategy);
	}

	void performDFS(GraphContext *g);
	void updateEdge(GraphContext *g, int u ,int v);
	void updateEdge(GraphContext *g, std::list<std::pair<int,int>> &edgeList);
	void reset();
	void displayTopologicalOrder();
	void displayReverseTopologicalOrder();
	void displayPostorder();
	void displayReversePostorder();
	void createMinI(GraphContext *g);
	void createMaxO(GraphContext *g);
	void displayMinI();
	void displayMaxO();
};
#endif