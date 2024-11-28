#ifndef ADJACENCYLISTGRAPH
#define ADJACENCYLISTGRAPH
#include"graph.h"
/*
Representation of the graph in adjacency list format inherting the necessary member 
functions and data members of the Graph class
*/
class AdjacencyListGraph : public Graph{
public:
    AdjacencyListGraph(int m):Graph(m) {}
    void addEdge(int x,int y) override;
    void addNode(int x) override;
    void printGOutgoing() const override;
    void printGIncoming() const override;
    const std::list<int>& getNeighboursIncoming(int x) const override;
    const std::list<int>& getNeighboursOutgoing(int x) const override;
    int getNumberOfVertices() const override;
    int getNumberOfEdges() const override;
};
#endif
