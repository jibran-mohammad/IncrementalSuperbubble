#ifndef SUPERBUBBLE_H1
#define SUPERBUBBLE_H1
#include<cstdio>
#include<stack>
#include<list>
#include<vector>
#include<climits>
#include<cassert>
#include<map>
#include<algorithm>
#include <utility>
#include <unordered_map>
#include <cstdint>
#include <limits>
#include"../graph/contextGraph.h"
#include"../DFS/contextdfs.h"


#define MAX(X,Y) ((X)>(Y)? (X):(Y))
#define MIN(X,Y) ((X)<(Y)? (X):(Y))
#define pii std::pair<int,int>
#define ppii std::pair<int,pii>

using SBID = uint32_t;

static constexpr SBID INVALID_SBD = std::numeric_limits<SBID>::max();
struct SuperbubbleNode {
    SBID id;
    SBID parent;                
    std::vector<SBID> children;
    int level;

    SuperbubbleNode(SBID _id)
        : id(_id), parent(INVALID_SBD), level(1) {}
};

class SuperbubbleCalculationH1 {
public:
    using SBID = uint32_t; 
    std::unordered_map<SBID, pii > SB; 
    std::unordered_map<SBID, SuperbubbleNode> SBTree;
    std::vector<std::vector<pii>> innerSB;
    std::vector<std::vector<int>> nodeList;
    std::vector<SBID> sbAtStart;
    std::vector<SBID> sbOfVertex;
private:
    SBID nextSBID;
public:
    SuperbubbleCalculationH1() : nextSBID(0) {}    
    void computeSB(DFSCONTEXT *sdfs,GraphContext *g);
    void printSuperbubbles() ;
    int numberOfSuperbubbles();
    void insertE(DFSCONTEXT *sdfs,GraphContext *g,int u, int v);
    void recomputeInside(DFSCONTEXT *sdfs,GraphContext *g,int s, int t);
private:
    SBID findLCA(SBID u, SBID v); 
    void clearSubtree(SBID rootID);
    void collectSubtree(SBID rootID, std::vector<SBID> &subtree);
       

};
#endif