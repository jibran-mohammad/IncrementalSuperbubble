#include <iostream>
#include <vector>
#include <list>
#include <random>
#include <fstream>
#include <algorithm>
#include <iostream>
#include "../graph/contextGraph.h"
#include "../DFS/contextdfs.h"
#include "../findSuperbubble/superbubble.h"
#include "../DFS/simpleDFS.h"
#include "../graph/adjacencyListGraph.h"

struct SuperbubbleSeedInfo {
    int seed;
    int edges;
    int numEdges;
};

class BarabasiAlbertDAG {
public:
    BarabasiAlbertDAG(GraphContext* g,int numNodes, int numEdges, int sd) 
        : g(g),numNodes(numNodes), numEdges(numEdges), seed(sd) {
            degrees.resize(numNodes, 0);
    }

    void generateDAG() {
        // Start with a small connected subgraph, avoiding nodes 0 and 1
        for (int i = 3; i < numEdges + 2; ++i) {
            g->insertEdge(2, i);
            degrees[2]++;
            degrees[i]++;
        }

        std::random_device rd;
        std::mt19937 gen(seed);

        // Add each new node from node 2 onwards
        for (int newNode = numEdges + 2; newNode < numNodes; ++newNode) {
            std::vector<int> targets = selectTargets(gen, newNode);

            for (int target : targets) {
                g->insertEdge(target, newNode);
                degrees[target]++;
                degrees[newNode]++;
            }
        }
    }


private:
    GraphContext* g; 
    int numNodes;
    int numEdges;
    int seed;
    std::vector<int> degrees;

    std::vector<int> selectTargets(std::mt19937& gen, int newNode) {
        std::vector<int> targets;
        std::discrete_distribution<> distribution(degrees.begin() + 2, degrees.begin() + newNode);

        while (targets.size() < numEdges) {
            int target = distribution(gen) + 2;  // Adjust for the offset, skipping nodes 0 and 1
            if (target != newNode && 
                std::find(targets.begin(), targets.end(), target) == targets.end()) {
                targets.push_back(target);
            }
        }
        return targets;
    }
};

int main() {
    long int trials = 150;
    int nodes = 10000;
    std::vector<SuperbubbleSeedInfo> superbubbleSeeds;
    int maxSuperbubbles = 0;
    int maxSuperbubblesSeed = -1;
    int maxSuperbubblesEdges = 0;
    int maxSuperbubblesNumEdges = 0;

    for (int t = 0; t < trials; ++t) {
        std::cout << "Trial: " << t << std::endl;
        
        // Generate a random number of edges within a moderate to dense range
        int edges = (rand() % static_cast<int>(0.2 * 49995000)) + static_cast<int>(0.3 * 49995000);
        // Vary numEdges based on edges and nodes, from sparse to more connected DAGs
        for (int numEdges = std::max(1, edges / nodes); numEdges <= std::min(nodes - 1, edges / nodes * 3); numEdges += std::max(1, edges / (15 * nodes))) {
            int seed = t;
            std::cout << "numEdges: " << numEdges << std::endl;
            AdjacencyListGraph adl(edges);
            GraphContext g(&adl);
            g.insertEdge(0,1);
            // Pass numEdges to the DAG generator
            BarabasiAlbertDAG dag(&g, nodes, numEdges, seed);
            dag.generateDAG();
            
            // Perform DFS and detect superbubbles
            SimpleDFS simdfs(&g);
            DFSCONTEXT d(&simdfs);
            SuperbubbleCalculation superbubble;
            
            d.performDFS(&g);
            d.createMinI(&g); 
            d.createMaxO(&g);
            superbubble.computeSB(&d, &g);

            // Calculate the number of superbubbles
            int superbubbleCount = superbubble.numberOfSuperbubbles();
            
            // Filter out graphs containing edges (0, 1)
            if (superbubbleCount > 0) {
                for (const auto& edge : superbubble.SB) {
                    if (edge.first == 0 && edge.second == 1) {
                        superbubbleCount = 0;
                        break;
                    }
                }
            }

            // Store seeds if superbubbles were found
            if (superbubbleCount > 0) {
                superbubbleSeeds.push_back({seed, edges, numEdges});
                
                // Update maximum superbubble count details
                if (superbubbleCount > maxSuperbubbles) {
                    maxSuperbubbles = superbubbleCount;
                    maxSuperbubblesSeed = seed;
                    maxSuperbubblesEdges = edges;
                    maxSuperbubblesNumEdges = numEdges;
                }
            }
        }
    }

    // Output results
    double superbubbleProbability = static_cast<double>(superbubbleSeeds.size()) / trials;
    std::cout << "Estimated superbubble probability: " << superbubbleProbability << std::endl;
    std::cout << "Seed with most superbubbles: " << maxSuperbubblesSeed 
              << " (superbubble count: " << maxSuperbubbles 
              << ", edges: " << maxSuperbubblesEdges 
              << ", numEdges: " << maxSuperbubblesNumEdges << ")" << std::endl;

    std::cout << "Seeds with superbubbles:" << std::endl;
    for (const auto& info : superbubbleSeeds) {
        std::cout << "(Seed: " << info.seed << ", Edges: " << info.edges << ", numEdges: " << info.numEdges << ") ";
    }
    std::cout << std::endl;

    return 0;
}