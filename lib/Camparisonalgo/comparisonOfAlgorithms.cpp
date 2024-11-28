#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>

// Function to read data from file and calculate total time
double calculateTotalTime(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return -1.0;
    }
    
    // Read the first line (number of edges and nodes)
    int numEdges, numNodes;
    file >> numEdges >> numNodes;
    
    // Variables to store edge data and total time
    int edgeNumber;
    double time;
    std::string caseType;
    double totalTime = 0.0;
    
    // Read the rest of the lines
    while (file >> edgeNumber >> time >> caseType) {
        totalTime += time; // Sum up time for each entry
    }
    
    file.close();
    return totalTime;
}

// Function to calculate the percentage speedup
double calculatePercentageSpeedup(double timeA, double timeB) {
    return ((timeA - timeB) / timeA) * 100.0;
}

// Function to calculate the speedup factor
double calculateSpeedupFactor(double timeA, double timeB) {
    return timeA / timeB;
}

int main() {
    // File paths for three algorithms
    std::string fileA = "/home/dell/Research/IncSuperbubble-master/output/Yeast/SuperbubbleSimpleTimeOutput.txt";
    std::string fileB = "/home/dell/Research/IncSuperbubble-master/output/Yeast/SuperbubbleincrementalDFSOutput.txt";
    std::string fileC = "/home/dell/Research/IncSuperbubble-master/output/Yeast/SuperbubbleincrementalDFSAVLTreeOutput.txt";
    
    // Calculate total times for each algorithm
    double totalTimeA = calculateTotalTime(fileA);
    double totalTimeB = calculateTotalTime(fileB);
    double totalTimeC = calculateTotalTime(fileC);

    if (totalTimeA < 0 || totalTimeB < 0 || totalTimeC < 0) {
        std::cerr << "Error in file processing." << std::endl;
        return 1;
    }

    // Display total times
    std::cout << "Total Time for Algorithm A: " << totalTimeA << " seconds\n";
    std::cout << "Total Time for Algorithm B: " << totalTimeB << " seconds\n";
    std::cout << "Total Time for Algorithm C: " << totalTimeC << " seconds\n";

    // Calculate and display speedup comparisons
    std::cout << "\nComparing Algorithm C and A:\n";
    double percentageSpeedupCA = calculatePercentageSpeedup(totalTimeC, totalTimeA);
    double speedupFactorCA = calculateSpeedupFactor(totalTimeC, totalTimeA);
    std::cout << "Percentage Speedup (C to A): " << percentageSpeedupCA << "%\n";
    std::cout << "Multiplicative Speedup Factor (A to B): " << speedupFactorCA << "x\n";

    std::cout << "\nComparing Algorithm C and B:\n";
    double percentageSpeedupCB = calculatePercentageSpeedup(totalTimeC, totalTimeB);
    double speedupFactorCB = calculateSpeedupFactor(totalTimeC, totalTimeB);
    std::cout << "Percentage Speedup (C to B): " << percentageSpeedupCB << "%\n";
    std::cout << "Multiplicative Speedup Factor (A to C): " << speedupFactorCB << "x\n";
    return 0;
}
