#ifndef GRAPH_H
#define GRAPH_H
#include<vector>
#include<list>
#include<map>
#include<algorithm>
#include<iterator>

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
	// for a vertex whose only parent is the artificial source 0: its place in
	// adjLOutgoing[0], so that the edge 0 -> y is removed in O(1) when y gets
	// a real parent (a std::list iterator stays valid until its element is
	// erased); list.remove(y) scanned all of 0's children
	std::vector<std::list<int>::iterator> sourcePos;

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

	//sizeV is the size of vertices and sizeE is the size of edges: every edge
	//in the graph, real and artificial (0->x, y->1, 0->1), counted by addEdge.

    int sizeE, sizeV;
    bool newu, newv;

    // Set by addEdge when the inserted edge u->v gave v its first real
    // parent, so the artificial edge 0->v was deleted (v was a child of the
    // artificial source in the DFS tree). Reset at the start of every addEdge.
    bool artificialParentRemoved;
    std::vector<int> incoming;
    std::vector<int> outgoing;
    // m (the number of edges of the input) is not the start value of sizeE:
    // addEdge counts every edge it adds, also the base graph's. (sizeE used
    // to start at m and count only the artificial edges; with the base graph
    // as m, the real edges of all later updates were never counted.)
    Graph(int m):sizeE(0),sizeV(0),newu(false),newv(false),artificialParentRemoved(false) {}

public:
	virtual void addEdge(int x,int y) = 0;
	virtual void addNode(int x) = 0;
	virtual void printGOutgoing() const = 0;
	virtual void printGIncoming() const = 0;
    virtual const std::list<int>& getNeighboursIncoming(int x) const =0;
    virtual const std::list<int>& getNeighboursOutgoing(int x) const =0;
    virtual int getNumberOfVertices() const = 0;
    virtual int getNumberOfEdges() const = 0;
    virtual int getInDegree(int x) const = 0;
	virtual int getOutDegree(int x) const = 0;
	virtual int getDegree(int x) const = 0;
    virtual ~Graph() = default;
};



#endif