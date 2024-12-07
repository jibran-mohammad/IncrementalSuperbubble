#include <iostream>
#include <fstream>
#include <vector>
#include <stack>
#include <string>
#include <unordered_set>
#include "../graph/adjacencyListGraph.h"
#include "../graph/contextGraph.h"

class ExtractDAGS {
public:
    ExtractDAGS(GraphContext &g) : visited(g.displayNumberVertices(), false), discoveryTime(g.displayNumberVertices(), -1),
                                   lowLink(g.displayNumberVertices(), -1), inStack(g.displayNumberVertices(), false), time(0) {}

    void findDAGs(const std::string& filename, GraphContext &g) {
        // Step 1: Find SCCs using Tarjan's algorithm
        std::vector<std::vector<int>> SCCs = findSCCs(g);

        int count = 1;
        int maxEdgeCount = 0;
        std::vector<std::pair<int, int>> largestDAGEdges;

        // Step 2: Extract and store each DAG separately, tracking the largest by edge count
        for (const auto& scc : SCCs) {
            // Only proceed if SCC has more than one node or at least one edge
            std::vector<std::pair<int, int>> dagEdges = getDAGEdges(scc, g);
            if (scc.size() > 1 || !dagEdges.empty()) {
                std::string outputFilename = filename + "_DAG" + std::to_string(count++) + ".el";
                std::ofstream outFile(outputFilename);
                if (!outFile.is_open()) {
                    std::cerr << "Error opening file for DAG output!" << std::endl;
                    continue;
                }

                outFile << dagEdges.size() << "\n";  // Number of edges
                for (const auto& edge : dagEdges) {
                    outFile << g.strategy->mapping[edge.first] << " " << g.strategy->mapping[edge.second] << "\n";
                }
                outFile.close();

                // Update largest DAG if this one has more edges
                if (dagEdges.size() > maxEdgeCount) {
                    maxEdgeCount = dagEdges.size();
                    largestDAGEdges = std::move(dagEdges);  // Move current edges to save memory
                }
            }
        }

        // Step 3: Save the largest DAG to a separate file
        if (!largestDAGEdges.empty()) {
            std::string largestFilename = filename + "_largest_DAG.el";
            std::ofstream largestOutFile(largestFilename);
            if (!largestOutFile.is_open()) {
                std::cerr << "Error opening file for largest DAG output!" << std::endl;
                return;
            }

            largestOutFile << maxEdgeCount << "\n";  // Number of edges in the largest DAG
            for (const auto& edge : largestDAGEdges) {
                largestOutFile << g.strategy->mapping[edge.first] << " " << g.strategy->mapping[edge.second] << "\n";
            }
            largestOutFile.close();
            std::cout << "Largest DAG with " << maxEdgeCount << " edges saved to " << largestFilename << std::endl;
        }
    }

private:
    std::vector<bool> visited;
    std::vector<int> discoveryTime, lowLink;
    std::vector<bool> inStack;
    std::stack<int> stack;
    int time;

    // Tarjan's algorithm to find SCCs
    std::vector<std::vector<int>> findSCCs(GraphContext &g) {
        std::vector<std::vector<int>> SCCs;
        for (int i = 0; i < g.displayNumberVertices(); ++i) {
            if (discoveryTime[i] == -1) {
                tarjanDFS(i, SCCs, g);
            }
        }
        return SCCs;
    }

    // DFS in Tarjan's algorithm
    void tarjanDFS(int u, std::vector<std::vector<int>>& SCCs, GraphContext &g) {
        discoveryTime[u] = lowLink[u] = time++;
        stack.push(u);
        inStack[u] = true;
        for (int v : g.displayOutgoingNeighbours(u)) {
            if (discoveryTime[v] == -1) {
                tarjanDFS(v, SCCs, g);
                lowLink[u] = std::min(lowLink[u], lowLink[v]);
            } else if (inStack[v]) {
                lowLink[u] = std::min(lowLink[u], discoveryTime[v]);
            }
        }

        // If u is a root node, pop all nodes in this SCC
        if (lowLink[u] == discoveryTime[u]) {
            std::vector<int> scc;
            while (true) {
                int v = stack.top();
                stack.pop();
                inStack[v] = false;
                scc.push_back(v);
                if (v == u) break;
            }
            SCCs.push_back(scc);
        }
    }

    // Get edges within SCC subgraph, ensuring acyclicity
    std::vector<std::pair<int, int>> getDAGEdges(const std::vector<int>& scc, GraphContext &g) {
        std::vector<std::pair<int, int>> dagEdges;
        std::unordered_set<int> nodeSet(scc.begin(), scc.end());
        
        for (int u : scc) {
            for (int v : g.displayOutgoingNeighbours(u)) {
                // Only include edges that do not form cycles within this SCC
                if (nodeSet.count(v) && u < v)  // Ensure acyclicity by restricting to u < v
                    dagEdges.push_back({u, v});
            }
        }
        return dagEdges;
    }
};

int main() {
    int m, x, y;
    std::cin >> m;
    AdjacencyListGraph adl(m);
    GraphContext g(&adl);
    g.insertEdge(0, 1);
    for (int i = 1; i <= m; i++) {
        std::cin >> x >> y;
        std::cout << x << " " << y << std::endl;
        g.insertEdge(x, y);
    }
    std::cout << "graph successfully created" << std::endl;

    // Initialize the graph with in-memory adjacency list
    ExtractDAGS graph(g);

    // Extract and store DAGs
    std::string filename = "/home/dell/Research/IncSuperbubble/output/webNotreDame/DAGS/web";  // Base filename
    graph.findDAGs(filename, g);

    std::cout << "DAGs extracted and saved to files." << std::endl;
    return 0;
}


