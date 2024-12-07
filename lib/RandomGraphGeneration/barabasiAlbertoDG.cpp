#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <string>
#include <algorithm>

class BarabasiAlbertDirected {
public:
    BarabasiAlbertDirected(int numNodes, int initialNodes, int edgesPerNode, int seed)
        : numNodes(numNodes), initialNodes(initialNodes), edgesPerNode(edgesPerNode), seed(seed) {
        degrees.resize(numNodes, 0);
    }

    void generateGraph() {
        std::random_device rd;
        std::mt19937 gen(seed);

        // Step 1: Initialize with a small fully connected graph
        for (int i = 0; i < initialNodes; ++i) {
            for (int j = 0; j < initialNodes; ++j) {
                if (i != j) {
                    std::cout << i <<" " << j << std::endl;
                    addEdge(i, j); // Fully connected
                }
            }
        }
        std::cout << "outside first for loop" << std::endl;
        // Step 2: Add new nodes with preferential attachment
        for (int newNode = initialNodes; newNode < numNodes; ++newNode) {
            std::vector<int> targets = selectTargets(gen, newNode);
            std::cout << "after select targets" << std::endl;
            for (int target : targets) {
                if (randomDirection(gen)) {
                    std::cout << newNode << " " << target << std::endl;
                    addEdge(newNode, target);
                } else {
                    std::cout << newNode << " " << target << std::endl;
                    addEdge(target, newNode);
                }
            }
        }
    }

    void saveGraphToFile(std::string& filename) const {
        std::string baseDir = "/home/dell/Research/IncSuperbubble/inputGraphs/BADirectedGraph/";
        filename = baseDir + filename;
        std::ofstream outFile(filename);
        if (!outFile.is_open()) {
            std::cerr << "Error: Unable to open file " << filename << std::endl;
            return;
        }
        outFile << edges.size() << "\n";
        for (const auto& edge : edges) {
            outFile << edge.first << " " << edge.second << "\n";
        }
        outFile.close();
        std::cout << "Graph saved to " << filename << std::endl;
    }

private:
    int numNodes;
    int initialNodes;
    int edgesPerNode;
    int seed;
    std::vector<int> degrees;
    std::vector<std::pair<int, int>> edges;

    void addEdge(int u, int v) {
        edges.emplace_back(u, v);
        degrees[u]++;
        degrees[v]++;
    }

    std::vector<int> selectTargets(std::mt19937& gen, int newNode) {
        std::vector<int> targets;
        std::discrete_distribution<> distribution(degrees.begin(), degrees.begin() + newNode);

        while (targets.size() < edgesPerNode) {
            int target = distribution(gen);
            if (target != newNode && std::find(targets.begin(), targets.end(), target) == targets.end()) {
                targets.push_back(target);
            }
        }
        return targets;
    }

    bool randomDirection(std::mt19937& gen) {
        std::uniform_int_distribution<> dist(0, 1);
        return dist(gen) == 1;
    }
};

int main() {
    const int numNodes = 10000;       // Number of nodes
    const int initialNodes = 5;       // Initial fully connected nodes
    const std::vector<int> edgeCounts = {5, 10, 15, 20}; // Edges per node to test
    const int numGraphsPerEdgeCount = 5; // Number of graphs per edge count

    for (int edgesPerNode : edgeCounts) {
        for (int seed = 0; seed < numGraphsPerEdgeCount; ++seed) {
            std::cout << "Generating graph with " << edgesPerNode 
                      << " edges per node and seed " << seed << "..." << std::endl;

            BarabasiAlbertDirected graph(numNodes, initialNodes, edgesPerNode, seed);
            graph.generateGraph();

            // Save the graph to a file
            std::string filename = "graph_" + std::to_string(edgesPerNode) + "edges_seed" 
                                    + std::to_string(seed) + ".txt";
            graph.saveGraphToFile(filename);
        }
    }

    std::cout << "All graphs generated and saved." << std::endl;
    return 0;
}
