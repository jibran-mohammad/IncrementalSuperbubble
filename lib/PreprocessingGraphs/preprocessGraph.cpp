#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

void processGraphFile(const std::string& inputFile, const std::string& outputFile) {
    std::ifstream inFile(inputFile);
    std::ofstream outFile(outputFile);

    if (!inFile.is_open()) {
        std::cerr << "Error: Unable to open input file: " << inputFile << std::endl;
        return;
    }
    
    if (!outFile.is_open()) {
        std::cerr << "Error: Unable to open output file: " << outputFile << std::endl;
        return;
    }

    std::string line;
    int numEdges = 0;

    // Read the first line (which is the number of edges)
    if (std::getline(inFile, line)) {
        //convert string into integer
        numEdges = std::stoi(line);
    }

    std::vector<std::pair<int, int>> edges;
    int validEdges = 0;

    // Read the remaining lines, process edges
    while (std::getline(inFile, line)) {
        std::istringstream iss(line);
        int u, v;

        // Check if the line contains valid integers
        if (iss >> u >> v) {
            // Skip lines where either u or v is 0 or 1
            if (u != 0 && v != 0 && u != 1 && v != 1) {
                edges.push_back({u, v});
                validEdges++;
            }
        }
    }

    // Output the number of valid edges
    outFile << validEdges << std::endl;

    // Output the filtered edges
    for (const auto& edge : edges) {
        outFile << edge.first << " " << edge.second << std::endl;
    }

    // Close the files
    inFile.close();
    outFile.close();

    std::cout << "Graph processing complete. Valid edges written to " << outputFile << std::endl;
}

int main() {
    std::string inputFile = "/home/dell/Research/IncSuperbubble/inputGraphs/webNotreDame/web-NotreDame.txt";  // Replace with your input file path
    std::string outputFile = "/home/dell/Research/IncSuperbubble/inputGraphs/webNotreDame/web-NotreDameProcessed.txt"; // Replace with your output file path

    processGraphFile(inputFile, outputFile);

    return 0;
}
