#ifndef CONTEXTGRAPH_H
#define CONTEXTGRAPH_H
#include"graph.h"
#include"adjacencyListGraph.h"
class GraphContext {
public:	
    Graph *strategy;

public:
	GraphContext(Graph *strategy): strategy(strategy){}
    
    void setStrategy(Graph *strategy){
		this->strategy = strategy;
	}
	void insertEdge(int x, int y);
	void insertNode(int x);
	void displayOutgoingGraph();
	void displayIncomingGraph();
	const std::list<int>& displayIncomingNeighbours(int);
	const std::list<int>& displayOutgoingNeighbours(int);
	int displayNumberEdges();
	int displayNumberVertices();
};	
#endif