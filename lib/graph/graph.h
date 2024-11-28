#ifndef GRAPH_H
#define GRAPH_H
#include<vector>
#include<list>
#include<map>
#include<algorithm>

/*
Graph interface defining necessary member functions for the graph and 
necessary data members
*/
class Graph{

public:
	//adjacency List for outgoing edges.

	std::vector<std::list<int> > adjLOutgoing;

	//adjacency List for incoming edges.

	std::vector<std::list<int> > adjLIncoming;

	/*
	array storing vertices of graph, indices will represent a 
	vertices from 0 to n-1, where n is number of vertices till now
	*/

	std::vector<int> mapping;

	/*
	key in the map is the original vertex, and the value is the 
    vertex we represent it with in mapping.
    */	

	std::map<int, int> reverseMapping;

	//sizeV is the size of vertices and sizeE is the size of edges.

    int sizeE, sizeV;
    std::vector<int> incoming;
    std::vector<int> outgoing;
    Graph(int m):sizeE(m),sizeV(0) {}

public:
	virtual void addEdge(int x,int y) = 0;
	virtual void addNode(int x) = 0;
	virtual void printGOutgoing() const = 0;
	virtual void printGIncoming() const = 0;
    virtual const std::list<int>& getNeighboursIncoming(int x) const =0;
    virtual const std::list<int>& getNeighboursOutgoing(int x) const =0;
    virtual int getNumberOfVertices() const = 0;
    virtual int getNumberOfEdges() const = 0;
    virtual ~Graph() = default;
};



#endif