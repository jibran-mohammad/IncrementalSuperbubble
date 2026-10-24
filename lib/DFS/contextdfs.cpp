#include "contextdfs.h"

    void DFSCONTEXT::performDFS(GraphContext *g) {
	    this->dfsStrategy->dfs(g);
    }

    void DFSCONTEXT::updateEdge(GraphContext *g, int u ,int v) {
	    this->dfsStrategy->insert(g,u,v);
    }

    void DFSCONTEXT::updateEdgeBatch(GraphContext *g, const std::vector<std::pair<int,int>> &edges) {
        this->dfsStrategy->insertBatch(g, edges);
    }

    void DFSCONTEXT::reset() {
	    this->dfsStrategy->reset();
    }

    void DFSCONTEXT::displayTopologicalOrder() {
	    this->dfsStrategy->printTopologicalOrder();
    }

    void DFSCONTEXT::displayReverseTopologicalOrder() {
	    this->dfsStrategy->printReverseTopologicalOrder();
    }

    void DFSCONTEXT::createMinI(GraphContext *g) {
	    this->datastructureStrategy->calculateMinI(g);
    }

    void DFSCONTEXT::createMaxO(GraphContext *g) {
	    this->datastructureStrategy->calculateMaxO(g);
    }

    void DFSCONTEXT::displayMinI() {
        this->datastructureStrategy->printMinI();
    }

    void DFSCONTEXT::displayMaxO() {
	    this->datastructureStrategy->printMaxO();
    }