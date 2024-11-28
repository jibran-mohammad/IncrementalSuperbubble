#include <iostream>
#include <vector>
#include <cstdlib>  // for rand() and srand()
#include <ctime>    // for time()
#include <set>
#include <algorithm> // for shuffle
#include <fstream>   // for file I/O

using namespace std;

// Function to generate a random DAG with 'V' vertices and 'E' edges, and save to file
void generateRandomDAG(int V, int E, double probability, const string& filename) {
    vector<pair<int, int>> edges; // Vector to store the edges
    set<pair<int, int>> uniqueEdges; // To ensure no duplicate edges

    srand(time(0)); // Initialize random seed

    // Step 1: Generate edges randomly between vertices
    for (int u = 0; u < V - 1; ++u) {
        for (int v = u + 1; v < V; ++v) {
            // Step 2: Only add an edge with a certain probability
            if ((rand() / (double)RAND_MAX) < probability) {
                std::cout << u << " " << v << std::endl;
                edges.push_back({u, v});
            }
        }
    }

    // Step 3: Shuffle the edges randomly to avoid sequentiality
    random_shuffle(edges.begin(), edges.end());

    // Step 4: Select the first 'E' edges
    if (E > edges.size()) {
        cout << "Too many edges requested for this number of vertices!" << endl;
        return;
    }

    // Step 5: Insert edges into the graph
    for (int i = 0; i < E; ++i) {
        uniqueEdges.insert(edges[i]);
    }

    // Step 6: Open the file for writing the graph
    ofstream outfile(filename);
    if (!outfile.is_open()) {
        cerr << "Error opening file!" << endl;
        return;
    }

    // Step 7: Write the number of vertices and edges in the first line
    outfile << E << endl;

    // Step 8: Write each edge to the file (u v)
    for (const auto& edge : uniqueEdges) {
        outfile << edge.first << " " << edge.second << endl;
    }

    // Close the file
    outfile.close();

    cout << "Graph with " << V << " vertices and " << E << " edges saved to " << filename << endl;
}

int main() {
    int V = 10000;  // Number of vertices (for example)
    int E = 1000000;  // Number of edges (over 1 million edges)
    double probability = 0.4;  // Probability of edge generation

    // Filename where the graph will be saved
    string filename = "/home/dell/Research/IncSuperbubble/inputGraphs/RandomAcyclicGraphs/randomgraph1.txt";

    // Generate random DAG with V vertices and E edges and save to file
    generateRandomDAG(V, E, probability, filename);

    return 0;
}
