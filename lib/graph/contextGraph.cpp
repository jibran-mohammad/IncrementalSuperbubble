#include"contextGraph.h"

    void GraphContext::insertEdge(int x, int y) {
    	this->strategy->addEdge(x,y);
    }

    void GraphContext::insertNode(int x) {
    	this->strategy->addNode(x);
    }

    void GraphContext::displayOutgoingGraph() {
    	this->strategy->printGOutgoing();
    }

    void GraphContext::displayIncomingGraph() {
    	this->strategy->printGIncoming();
    }

    const std::list<int>& GraphContext::displayIncomingNeighbours(int x) {
    	this->strategy->getNeighboursIncoming(x);
    }

    const std::list<int>& GraphContext::displayOutgoingNeighbours(int x) {
    	this->strategy->getNeighboursOutgoing(x);
    }

    int GraphContext::displayNumberEdges() {
    	this->strategy->getNumberOfEdges();
    }

    int GraphContext::displayNumberVertices() {
    	this->strategy->getNumberOfVertices();
    }
