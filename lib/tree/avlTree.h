#ifndef AVLTREE_
#define AVLTREE_
#include<map>
#include<iostream>
/*
  The NORMAL balanced binary search tree (std::map, a red-black tree) used by
  incNormalBST: the key is a vertex's position in the topological order
  (reverseTopologicalOrder), the value the vertex (index). findMin / findMax
  give the vertex with the smallest / largest position. When a vertex moves,
  its key must be deleted and added again with the new position (compare
  differentAVLTree.h, incBST's modified tree, which stores vertex ids and
  reads their current positions when comparing).
*/
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