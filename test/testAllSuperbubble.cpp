
#include <iostream>
#include <fstream>
#include <string>

void compareFiles(const std::string& file1, const std::string& file2,const std::string& file3, const std::string& diffFile) {
    std::ifstream infile1(file1);
    std::ifstream infile2(file2);
    std::ifstream infile3(file3);
    std::ofstream diffout(diffFile);

    // Check if files opened successfully
    if (!infile1.is_open()) {
        std::cerr << "Error: Unable to open " << file1 << std::endl;
        return;
    }
    if (!infile2.is_open()) {
        std::cerr << "Error: Unable to open " << file2 << std::endl;
        return;
    }
    if (!infile3.is_open()) {
        std::cerr << "Error: Unable to open " << file3 << std::endl;
        return;
    }
    if (!diffout.is_open()) {
        std::cerr << "Error: Unable to open " << diffFile << " for writing." << std::endl;
        return;
    }

    std::string line1, line2, line3;
    int lineNumber = 0;
    bool filesAreIdentical = true;

    while (true) {
        std::getline(infile1, line1);
        std::getline(infile2, line2);
        std::getline(infile3, line3);
        lineNumber++;

        // Check if all three files reached EOF
        if (infile1.eof() && infile2.eof() && infile3.eof()) {
            break; // all files are the same length and ended at the same time
        }

        // If one file reached EOF, and the others didn't
        if (infile1.eof() || infile2.eof() || infile3.eof()) {
            filesAreIdentical = false;
            diffout << "Difference at line " << lineNumber << ": ";
            if (infile1.eof()) {
                diffout << "End of " << file1 << ", content of " << file2 << ": " << line2 << std::endl;
            } 
            else if(infile2.eof()) {
                diffout << "End of " << file2 << ", content of " << file1 << ": " << line1 << std::endl;
            }
            else {
                diffout << "End of " << file3 << ", content of " << file3 << ": " << line3 << std::endl;
            }
            continue; // No need to compare lines further
        }

        // Compare the two lines
        if (line1 != line2 || line2 != line3 || line1 !=line3) {
            filesAreIdentical = false;
            diffout << "Difference at line " << lineNumber << ": \n";
            diffout << "File 1: " << line1 << "\n";
            diffout << "File 2: " << line2 << "\n";
            diffout << "File 3: " << line3 << "\n";
        }
    }

    if (filesAreIdentical) {
        std::cout << "Files are identical." << std::endl;
    } else {
        std::cout << "Files are different. Differences logged to " << diffFile << std::endl;
    }

    infile1.close();
    infile2.close();
    infile2.close();
    diffout.close();
}

int main() {
    std::string file1 = "/home/dell/Research/IncSuperbubble-master/output/Email/SuperbubblesimpleDFSOutput.txt"; // Change this to your file path
    std::string file2 = "/home/dell/Research/IncSuperbubble-master/output/Email/SuperbubbleincrementalDFSOutput.txt"; // Change this to your file path
    std::string file3 = "/home/dell/Research/IncSuperbubble-master/output/Email/SuperbubbleincrementalDFSAVLTreeOutput.txt";
    std::string diffFile = "/home/dell/Research/IncSuperbubble-master/output/Yeast/difference.txt"; // Change this to your desired output path

    compareFiles(file1, file2, file3, diffFile);

    return 0;
}