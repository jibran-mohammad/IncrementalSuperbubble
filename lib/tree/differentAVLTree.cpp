#include"differentAVLTree.h"

    // Comparison operator
    bool TopologicalCompare::operator()(const int &a, const int &b) const {
        // Compare the positions of a and b in the topological order
        return (*reversePostOrder)[a] < (*reversePostOrder)[b];
    }


    // Function to add a node to the binary tree
    void DifferentBalancedBinaryTree::addNode(int index, int reversePostOrder) {
        avlTree.insert(std::make_pair(index,reversePostOrder));
    }

    // Function to find the maximum element in the binary tree
    int DifferentBalancedBinaryTree::findMax() {
        if (avlTree.empty()){
            return 1;
        }
        else {    
           return avlTree.rbegin()->first;
        }
    }

    // Function to find the minimum element in the binary tree
    int DifferentBalancedBinaryTree::findMin() {
        if (avlTree.empty()){
           return 0; 
        }
        else {
           return avlTree.begin()->first;
        }
    }

    // Function to delete a node from the binary tree
    void DifferentBalancedBinaryTree::deleteNode(int index) {
        auto it = avlTree.find(index);
        if (it != avlTree.end()) {
            //std::cout << "inside deletes if"<<std::endl;
            avlTree.erase(it);
        }
        //else
            //std::cout << "Node not found" << std::endl;
    }

    // Function to print the binary tree
    void DifferentBalancedBinaryTree::printTree() {
        for (const auto& pair : avlTree) {
            std::cout << "index: " << pair.first << ", topologicalOrder: " << pair.second << std::endl;
        }
    }