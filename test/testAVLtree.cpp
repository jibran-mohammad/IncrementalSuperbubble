#include"../lib/tree/avlTree.h"

int main() {
    BalancedBinaryTree binaryTree;

    // Add nodes to the binary tree
    binaryTree.addNode(3, 10);
    binaryTree.addNode(5, 20);
    binaryTree.addNode(2, 30);
    binaryTree.addNode(7, 40);

    // Print the binary tree
    std::cout << "Binary Tree:" << std::endl;
    binaryTree.printTree();

    // Find maximum and minimum elements
    TreeNode maxNode = binaryTree.findMax();
    TreeNode minNode = binaryTree.findMin();
    std::cout << "Maximum: Topological Order - " << maxNode.reversetopologicalOrder << ", Index - " << maxNode.index << std::endl;
    std::cout << "Minimum: Topological Order - " << minNode.reversetopologicalOrder << ", Index - " << minNode.index << std::endl;

    // Delete a node from the binary tree
    binaryTree.deleteNode(5);
    std::cout << "Binary Tree after deleting node with topological order 5:" << std::endl;
    binaryTree.printTree();

    return 0;
}
