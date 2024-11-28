#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include "../graph/contextGraph.h"
#include "../DFS/contextdfs.h"
#include "../findSuperbubble/superbubble.h"
#include "../DFS/simpleDFS.h"
#include "../graph/adjacencyListGraph.h"
#include "../DFS/incrementalDFSAVLTree.h"
#include "../DFS/incrementalDFS.h"

struct SuperbubbleSeedInfo {
    int seed;
    int edges;
};

// Generate a random DAG based on a seed, excluding nodes 0 and 1
void generateRandomDAG(int nodes, int edges, int seed, GraphContext *g) {
    std::mt19937 rng(seed);
    std::vector<int> topoOrder(nodes - 2); // Exclude nodes 0 and 1

    for (int i = 2; i < nodes; ++i) {
        topoOrder[i - 2] = i;  // Shift by 2 to accommodate the exclusion of 0 and 1
    }

    // Shuffle the topological order
    std::shuffle(topoOrder.begin(), topoOrder.end(), rng);
    
    int edgeCount = 0;
    g->insertEdge(0,1);
    // Generate edges only between nodes in the topoOrder, avoiding nodes 0 and 1
    for (int i = 0; i < nodes - 2 && edgeCount < edges; ++i) {
        for (int j = i + 1; j < nodes - 2 && edgeCount < edges; ++j) {
            if (rng() % 2) {
                g->insertEdge(topoOrder[i], topoOrder[j]);
                edgeCount++;
            }
        }
    }
}

int main() {
    long int trials = 10000;
    int nodes = 1000;
    std::vector<SuperbubbleSeedInfo> superbubbleSeeds;
    int maxSuperbubbles = 0;
    int maxSuperbubblesSeed = -1;
    int maxSuperbubblesEdges = 0;
    for (int t = 0; t < trials; ++t) {
        std::cout << t << std::endl;
        //int edges = (rand() % 6) + 22; // 22 t0 28
        //int edges = (rand() % (int)(0.5 * 66)) + (0.9 * 66); // 17 to 22
        //int edges = (rand() % (59 - 33 + 1)) + 33;
        int edges = (rand() % static_cast<int>(0.3 * 499500)) + static_cast<int>(0.4 * 499500);
        int seed = t;
        AdjacencyListGraph adl(edges) ;
        GraphContext g(&adl);
        generateRandomDAG(nodes, edges, seed,&g);
        SimpleDFS simdfs(&g);
        DFSCONTEXT d(&simdfs);
        SuperbubbleCalculation superbubble;
        d.performDFS(&g);
        d.createMinI(&g);
        d.createMaxO(&g);
        superbubble.computeSB(&d,&g);

        int superbubbleCount = superbubble.numberOfSuperbubbles();
        
        if (superbubbleCount > 0) {
            for (const auto& edge : superbubble.SB) {
                if (edge.first == 0 && edge.second == 1) {
                  superbubbleCount -= superbubbleCount;
                  break;
                }
            }
        }

        if (superbubbleCount > 0) {
            superbubbleSeeds.push_back({seed, edges});
            if (superbubbleCount > maxSuperbubbles) {
                maxSuperbubbles = superbubbleCount;
                maxSuperbubblesSeed = seed;
            }
        }

    }

    double superbubbleProbability = static_cast<double>(superbubbleSeeds.size()) / trials;
    std::cout << "Estimated superbubble probability: " << superbubbleProbability << std::endl;
    std::cout << "Seed with most superbubbles: " << maxSuperbubblesSeed  << " (superbubble count: " << maxSuperbubbles << ", edges: " << maxSuperbubblesEdges << ")" << std::endl;

    std::cout << "Seeds with superbubbles:" << std::endl;
    for (const auto& info : superbubbleSeeds) {
        std::cout << "(" << info.seed << ", " << info.edges << ") ";
    }
    std::cout << std::endl;

    return 0;
}
