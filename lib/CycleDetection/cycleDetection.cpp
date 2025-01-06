#include"../graph/adjacencyListGraph.h"
#include"../graph/contextGraph.h"
#include<iostream>

bool dfs(int v, std::vector<bool>& visited, std::vector<bool>& recStack, GraphContext &g) {
        visited[v] = true;
        recStack[v] = true;

        // Recur for all adjacent vertices
        for (int neighbor : g.displayOutgoingNeighbours(v)) {
            // If neighbor is not visited, visit it
            if (!visited[neighbor] && dfs(neighbor, visited, recStack, g)) {
                return true;
            }
            // If neighbor is in the recursion stack, we found a back edge
            else if (recStack[neighbor]) {
                return true;
            }
        }

        // Remove vertex from recursion stack
        recStack[v] = false;
        return false;
}

bool isGraphDAG(GraphContext &g) {
        std::vector<bool> visited(g.displayNumberVertices(), false);
        std::vector<bool> recStack(g.displayNumberVertices(), false);

        // Perform DFS from every unvisited vertex
        for (int i = 0; i < g.displayNumberVertices() ; ++i) {
            if (!visited[i]) {
                if (dfs(i, visited, recStack, g)) {
                    return false;  // Cycle detected
                }
            }
        }
        return true;  // No cycle detected
}
    
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
    if (isGraphDAG(g)) {
        std::cout << "The graph is a DAG." << std::endl;
    } else {
        std::cout << "The graph is not a DAG." << std::endl;
    }
	return 0;


}