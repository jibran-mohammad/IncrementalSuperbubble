 #include"avlTree.h"


    // Function to add a node to the binary tree
    void BalancedBinaryTree::addNode(int reversetopOrder, int idx) {
        avlTree.insert(std::make_pair(reversetopOrder,TreeNode(reversetopOrder, idx)));
    }

    // Function to find the maximum element in the binary tree
    TreeNode BalancedBinaryTree::findMax() {
        TreeNode tree;
        if (avlTree.empty()){
           return tree;
        }
        else {    
           return avlTree.rbegin()->second;
        }
    }

    // Function to find the minimum element in the binary tree
    TreeNode BalancedBinaryTree::findMin() {
        TreeNode tree;
        if (avlTree.empty()){
           return tree; 
        }
        else {
           return avlTree.begin()->second;
        }
    }

    // Function to delete a node from the binary tree
    void BalancedBinaryTree::deleteNode(int reversetopOrder) {
        auto it = avlTree.find(reversetopOrder);
        if (it != avlTree.end())
            avlTree.erase(it);
    }

    // Function to print the binary tree
    void BalancedBinaryTree::printTree() {
        for (const auto& pair : avlTree) {
            std::cout << "Topological Order: " << pair.second.reversetopologicalOrder << ", Index: " << pair.second.index << std::endl;
        }
    }

