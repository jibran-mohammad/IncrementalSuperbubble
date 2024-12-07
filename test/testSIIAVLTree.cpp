#include <iostream>
#include <fstream>
#include <string>

void compareFiles(std::string& file1, std::string& file2, std::string& file3, std::string& file4, std::string& diffFile) {
    std::string baseDirInput1 = "/home/dell/Research/IncSuperbubble/output/HunderedRandomGraphs/SimpleDFS/";
    std::string baseDirInput2 = "/home/dell/Research/IncSuperbubble/output/HunderedRandomGraphs/IncrementalDFS/";
    std::string baseDirInput3 = "/home/dell/Research/IncSuperbubble/output/HunderedRandomGraphs/IncrementalDFSAVLTree/";
    std::string baseDirInput4 = "/home/dell/Research/IncSuperbubble/output/HunderedRandomGraphs/incrementalDFSH1/";
    std::string differenceOutput = "/home/dell/Research/IncSuperbubble/output/HunderedRandomGraphs/Difference/";
    file1 = baseDirInput1 + file1 ;
    file2 = baseDirInput2 + file2 ;
    file3 = baseDirInput3 + file3 ;
    file4 = baseDirInput4 + file4 ;
    diffFile = differenceOutput + diffFile; 

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

        // Check if all three files reached EOF
        if (infile1.eof() && infile2.eof() && infile3.eof() && infile4.eof()) {
            break; // all files are the same length and ended at the same time
        }

        // If one file reached EOF, and the others didn't
        if (infile1.eof() || infile2.eof() || infile3.eof() || infile4.eof()) {
            filesAreIdentical = false;
            diffout << "Difference at line " << lineNumber << ": ";
            if (infile1.eof()) {
                diffout << "End of " << file1 << ", content of " << file2 << ": " << line2 << std::endl;
            } 
            else if(infile2.eof()) {
                diffout << "End of " << file2 << ", content of " << file1 << ": " << line1 << std::endl;
            }
            else if(infile3.eof()){
                diffout << "End of " << file3 << ", content of " << file3 << ": " << line3 << std::endl;
            }
            else {
                diffout << "End of " << file4 << ", content of " << file4 << ": " << line4 << std::endl;
            }
            continue; // No need to compare lines further
        }

        // Compare the two lines
        if (line1 != line2 || line2 != line3 || line1 !=line3 || line1 != line4 || line2 != line4 || line3 != line4) {
            filesAreIdentical = false;
            diffout << "Difference at line " << lineNumber << ": \n";
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
    std::string filename1;
    std::string filename2;
    std::string filename3;
    std::string filename4;
    std::string diffFile;
    
    for(int i=1; i<=100; i++) {
        filename1 = "randomGraph" + std::to_string(i) + "simpleDFS.txt";
        filename2 = "randomGraph" + std::to_string(i) + "incrementalDFS.txt";
        filename3 = "randomGraph" + std::to_string(i) + "incrementalDFSAvlTree.txt";
        filename4 = "randomGraph" + std::to_string(i) + "incrementalDFSH1.txt";
        diffFile = "difference" + std::to_string(i) + ".txt";
        compareFiles(filename1, filename2, filename3,filename4, diffFile);
    }
    
    //filename1 = "randomgraph2processedsimpleDFS.txt";
    //filename2 = "randomgraph2processedincrementalDFS.txt";
    //filename3 = "randomgraph2processedincrementalDFSAvlTree.txt";
    //diffFile = "difference.txt";
    //compareFiles(filename1, filename2, filename3, diffFile);

    return 0;
}
