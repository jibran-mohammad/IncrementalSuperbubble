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
#include<ostream>
#include"../graph/contextGraph.h"
#include"../DFS/contextdfs.h"


#define MAX(X,Y) ((X)>(Y)? (X):(Y))
#define MIN(X,Y) ((X)<(Y)? (X):(Y))
#define pii std::pair<int,int>
#define ppii std::pair<int,pii>

class SuperbubbleCalculation {
public:
    std::list<pii> SB;
    // time of clearing SB in the last insertE / insertBatch (normal build):
    // the incremental programs add it to UpdateDS; static times the whole
    // insertE as UpdateDS, so there it is already included
    double clearTime = 0.0;
public:    
    void computeSB(DFSCONTEXT *sdfs,GraphContext *g);
    void printSuperbubbles(GraphContext *g) ;
    void writeSuperbubbles(GraphContext *g, std::ostream &out);
    int numberOfSuperbubbles();
    void insertE(DFSCONTEXT *sdfs,GraphContext *g,int u, int v);
    // the edges of a batch (input ids), applied at once by the DFS strategy
    void insertBatch(DFSCONTEXT *sdfs, GraphContext *g, const std::vector<std::pair<int,int>> &edges);
};
#endif	
