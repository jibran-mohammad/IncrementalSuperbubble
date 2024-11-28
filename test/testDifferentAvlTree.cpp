#include"../lib/tree/differentAVLTree.h"

int main() {
    std::vector<int> reversetopologicalOrder = {0,1,2,3,11,6,8,10,12,7,9,4,5,13,15,14,16};
    DifferentBalancedBinaryTree binaryTree(reversetopologicalOrder);  
    
    // Add nodes to the binary tree
    binaryTree.addNode(0, 0);
    binaryTree.addNode(1, 0);
    binaryTree.addNode(2, 0);
    binaryTree.addNode(3, 0);
    binaryTree.addNode(4, 0);
    binaryTree.addNode(5, 0);
    binaryTree.addNode(6, 0);
    binaryTree.addNode(7, 0);
    binaryTree.addNode(8, 0);
    binaryTree.addNode(9, 0);
    binaryTree.addNode(10, 0);
    binaryTree.addNode(11, 0);
    binaryTree.addNode(12, 0);
    binaryTree.addNode(13, 0);
    binaryTree.addNode(14, 0);
    binaryTree.addNode(15, 0);
    binaryTree.addNode(16, 0);
    // Print the binary tree
    std::cout << "Binary Tree:" << std::endl;
    binaryTree.printTree();

    // Find maximum and minimum elements
    int maxNode = binaryTree.findMax();
    int minNode = binaryTree.findMin();
    std::cout << "Maximum: node - " << maxNode << ", Index - " << std::endl;
    std::cout << "Minimum: node - " << minNode << ", Index - " << std::endl;
    // Delete a node from the binary tree
    binaryTree.deleteNode(5);
    std::cout << "after deleting 5" << std::endl;
    binaryTree.printTree();
    binaryTree.deleteNode(7);
    binaryTree.deleteNode(11);
    binaryTree.deleteNode(13);
    binaryTree.deleteNode(14);
    binaryTree.deleteNode(16);
    reversetopologicalOrder = {0,1,2,3,11,10,8,6,12,7,9,13,5,4,16,14,15};

    binaryTree.printTree();
    binaryTree.addNode(5,0);
    binaryTree.addNode(7,0);
    binaryTree.addNode(11,0);
    binaryTree.addNode(13,0);
    binaryTree.addNode(14,0);
    binaryTree.addNode(16,0);
    std::cout << "Binary Tree after changed topologicalOrder:" << std::endl;
    binaryTree.printTree();
    return 0;
}
