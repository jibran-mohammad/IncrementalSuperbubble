#include"../lib/graph/adjacencyListGraph.h"
#include"../lib/graph/contextGraph.h"
#include<iostream>


/*
test file creating the objects of adjacecnyListGraph class and calling 
appropriate functions for creating the graph and printing the graph 
as adjacency list with incoming and outgoing edges.
*/
int main()
{
	int m,x,y;
    std::cin >> m;
    AdjacencyListGraph adl(m) ;
    GraphContext g(&adl);
    g.insertEdge(0,1);
    for (int i = 1; i <= m; i++) {
           std::cin >> x >> y;
           g.insertEdge(x,y);
           if(x== 5966 && y== 2913) {
           	for(int elem:g.displayOutgoingNeighbours(16335)) {
		std::cout << elem << " " << std::endl;
	       }
           }
           /*
	       g.displayOutgoingGraph();
	       g.displayIncomingGraph();
	       std::cout << "Number of edges = " << g.displayNumberEdges() << std::endl;
	       std::cout << "Number of vertices = " << g.displayNumberVertices() << std::endl;
	     */
	}
	/*
	for(int elem:g.displayOutgoingNeighbours(16335)) {
		std::cout << elem << " " << std::endl;
	}
	*/
	/*
	g.displayOutgoingGraph();
	g.displayIncomingGraph();
	std::cout << "Number of edges = " << g.displayNumberEdges() << std::endl;
	std::cout << "Number of vertices = " << g.displayNumberVertices() << std::endl;
	
	*/
	return 0;


}