#include <iostream>
#include <vector>
#include <cstdlib>  // for rand() and srand()
#include <ctime>    // for time()
#include <set>
#include <algorithm> // for shuffle
#include <fstream>   // for file I/O
#include <string>    // for filenames

using namespace std;

// Function to generate a random DAG with 'V' vertices and 'E' edges, and save to file
void generateRandomDAG(int V, int E, double probability, const string& filename) {
    vector<pair<int, int>> edges; // Vector to store the edges
    set<pair<int, int>> uniqueEdges; // To ensure no duplicate edges

    // Step 1: Generate edges randomly between vertices
    for (int u = 2; u <V - 1; ++u) {
        for (int v = u + 1; v < V ; ++v) {
            // Step 2: Only add an edge with a certain probability
            if ((rand() / (double)RAND_MAX) < probability) {
                edges.push_back({u, v});
            }
        }
    }

    // Step 3: Shuffle the edges randomly to avoid sequentiality
    random_shuffle(edges.begin(), edges.end());

    // Step 4: Select the first 'E' edges
    if (E > edges.size()) {
        cout << "Requested " << E << " edges, but only " << edges.size() << " possible. Adjusting to " << edges.size() << " edges." << endl;
        E = edges.size();
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
}

// Main function to generate 100 random graphs with varying vertices and edges
int main() {
    srand(time(0));  // Initialize random seed

    int numGraphs = 100;  // Number of graphs to generate

    // Generate 100 random graphs with varying vertices and edges
    for (int i = 1; i <= numGraphs; ++i) {
        // Randomly choose the number of vertices between 10 and 15
        int V = rand() % 6 + 10;  // Random number between 10 and 15

        // Randomly choose the number of edges between 20 and 25
        int E = rand() % 6 + 20;  // Random number between 20 and 25

        // Set a random edge generation probability (optional adjustment)
        double probability = 0.5;  // Adjust this value if you want to control density

        // Create a unique filename for each graph
        string filename = "/home/dell/Research/IncSuperbubble-master/test/inputGraphs/hunderedRandomGraphs/randomGraph" + to_string(i) + ".txt";

        // Generate and save the random DAG
        generateRandomDAG(V, E, probability, filename);

        cout << "Graph " << i << " with " << V << " vertices and " << E << " edges saved to " << filename << endl;
    }

    return 0;
}
