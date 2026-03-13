#ifndef DFS
#define DFS
#include<vector>
#include<iostream>
#include<algorithm>
#include "../graph/contextGraph.h"
/*
 Base class for defining the functionality used by all the algorithms calculating the postorder, 
 reverpostorder, toporder, revertoporder.
*/
class IDFS{
public:	
	std::vector<int> visit;
    std::vector<int> topOrd;
    std::vector<int> postorder;
    std::vector<int> reversePostorder;
    std::vector<int> reverseTopologicalOrder;
    std::vector<int> oldReverseTopologicalOrder;
    std::vector<int> minI;
    std::vector<int> maxO;
    int n;
    int postorderNumber,cases;

public:
	
    virtual void dfs(GraphContext *g) = 0;
    virtual void printTopologicalOrder() const = 0;
    virtual void printReverseTopologicalOrder() const = 0;
    virtual void printPostOrder() const = 0;
    virtual void insert(GraphContext * g, int u , int v) = 0;
    virtual void printreversePostOrder() const = 0;
    virtual void reset() = 0;

    virtual ~IDFS() = default;
};
#endif