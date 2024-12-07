
#include <iostream>
#include <fstream>
#include <string>

void compareFiles(const std::string& file1, const std::string& file2, const std::string& file3, const std::string& file4, const std::string& diffFile) {
    std::ifstream infile1(file1);
    std::ifstream infile2(file2);
    std::ifstream infile3(file3);
    std::ifstream infile4(file4);
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
    if (!infile4.is_open()) {
        std::cerr << "Error: Unable to open " << file4 << std::endl;
        return;
    }
    if (!diffout.is_open()) {
        std::cerr << "Error: Unable to open " << diffFile << " for writing." << std::endl;
        return;
    }

    std::string line1, line2, line3, line4;
    int lineNumber = 0;
    bool filesAreIdentical = true;

    while (true) {
        std::getline(infile1, line1);
        std::getline(infile2, line2);
        std::getline(infile3, line3);
        std::getline(infile4, line4);
        lineNumber++;

        // Check if all four files reached EOF
        if (infile1.eof() && infile2.eof() && infile3.eof() && infile4.eof()) {
            break; // all files are the same length and ended at the same time
        }

        // If one file reached EOF, and the others didn't
        if (infile1.eof() || infile2.eof() || infile3.eof() || infile4.eof()) {
            filesAreIdentical = false;
            diffout << "Difference at line " << lineNumber << ": ";
            if (infile1.eof()) {
                diffout << "End of " << file1 << ", contents of other files continue." << std::endl;
            }
            if (infile2.eof()) {
                diffout << "End of " << file2 << ", contents of other files continue." << std::endl;
            }
            if (infile3.eof()) {
                diffout << "End of " << file3 << ", contents of other files continue." << std::endl;
            }
            if (infile4.eof()) {
                diffout << "End of " << file4 << ", contents of other files continue." << std::endl;
            }
            continue; // No need to compare lines further
        }

        // Compare the lines
        if (line1 != line2 || line2 != line3 || line3 != line4 || line1 != line4) {
            filesAreIdentical = false;
            diffout << "Difference at line " << lineNumber << ":\n";
            diffout << "File 1: " << line1 << "\n";
            diffout << "File 2: " << line2 << "\n";
            diffout << "File 3: " << line3 << "\n";
            diffout << "File 4: " << line4 << "\n";
        }
    }

    if (filesAreIdentical) {
        std::cout << "Files are identical." << std::endl;
    } else {
        std::cout << "Files are different. Differences logged to " << diffFile << std::endl;
    }

    infile1.close();
    infile2.close();
    infile3.close();
    infile4.close();
    diffout.close();
}

int main() {
    std::string file1 = "/home/dell/Research/IncSuperbubble/output/Yeast/superbubbles/SuperbubblesimpleDFSOutput.txt"; 
    std::string file2 = "/home/dell/Research/IncSuperbubble/output/Yeast/superbubbles/SuperbubbleincrementalDFSOutput.txt"; 
    std::string file3 = "/home/dell/Research/IncSuperbubble/output/Yeast/superbubbles/SuperbubbleincrementalDFSAVLTreeOutput.txt";
    std::string file4 = "/home/dell/Research/IncSuperbubble/output/Yeast/superbubbles/SuperbubbleincrementalH1output.txt"; 
    std::string diffFile = "/home/dell/Research/IncSuperbubble/output/Yeast/difference.txt"; 

    compareFiles(file1, file2, file3, file4, diffFile);

    return 0;
}
