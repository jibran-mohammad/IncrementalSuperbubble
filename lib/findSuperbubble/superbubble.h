#ifndef SUPERBUBBLE
#define SUPERBUBBLE
#include<cstdio>
#include<stack>
#include<list>
#include<vector>
#include<climits>
#include<cassert>
#include<map>
#include<algorithm>
#include"../graph/contextGraph.h"
#include"../DFS/contextdfs.h"


#define MAX(X,Y) ((X)>(Y)? (X):(Y))
#define MIN(X,Y) ((X)<(Y)? (X):(Y))
#define pii std::pair<int,int>
#define ppii std::pair<int,pii>

class SuperbubbleCalculation {
public:
    std::list<pii> SB;
public:    
    void computeSB(DFSCONTEXT *sdfs,GraphContext *g);
    void printSuperbubbles() ;
    int numberOfSuperbubbles();
    void insertE(DFSCONTEXT *sdfs,GraphContext *g,int u, int v);
    void insertBatch(DFSCONTEXT *sdfs, GraphContext *g);
};
#endif	
