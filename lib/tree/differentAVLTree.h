#ifndef DIFFERENT_AVL_TREE
#define DIFFERENT_AVL_TREE
#include<map>
#include<iostream>
#include<vector>

class TopologicalCompare {
public:
    std::vector<int> *reversePostOrder;
    TopologicalCompare() = default;

    // Constructor that takes the topological order vector
    TopologicalCompare(std::vector<int> &order) {
        reversePostOrder = &order;
    }
    bool operator()(const int &a, const int &b) const;
    
};
class DifferentBalancedBinaryTree {


    public:
        std::map<int,int,TopologicalCompare> avlTree;
        DifferentBalancedBinaryTree() = default;
        DifferentBalancedBinaryTree(std::vector<int> &order) : avlTree(TopologicalCompare(order)) {}
        void addNode(int index, int reversePostOrder);
        int findMax();
        int findMin();
        void deleteNode(int index);
        void printTree();
        void assignTopologicalOrder();
};
#endif