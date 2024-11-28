#ifndef AVLTREE_
#define AVLTREE_
#include<map>
#include<iostream>
struct TreeNode {
    int reversetopologicalOrder;
    int index;
    TreeNode() : reversetopologicalOrder(0), index(0) {}
    TreeNode(int reversetopOrder, int idx) : reversetopologicalOrder(reversetopOrder), index(idx) {}
};

class BalancedBinaryTree {
    public:
        std::map<int, TreeNode> avlTree;
        void addNode(int reversetopOrder, int idx);
        TreeNode findMax();
        TreeNode findMin();
        void deleteNode(int reversetopOrder);
        void printTree();
};
#endif