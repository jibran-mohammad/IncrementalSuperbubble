#include"../lib/findSuperbubble/superbubblesimpleDFS.h"

int main(){
	AdjacencyListGraph adl ;
	Graph *g = &adl;
	IDFS *dfs;
	g->createGInp();
	g->printGOutgoing();
	g->printGIncoming();
	SimpleDFS simdfs(g);
	dfs = &simdfs;
	dfs->dfs(g);
	dfs->printReversePostOrder(g);
	dfs->printTopologicalOrder(g);
	SuperbubbleCalculationSimpleDFS superc(dfs,g);
	superc.calculateMinI(dfs,g);
	superc.calculateMaxO(dfs,g);
	superc.printMaxO();
    superc.printMinI();
    superc.computeSB(dfs,g);
    
    superc.printSuperbubbles(dfs,g);
    std::cout << "The number of superbubbles are" << std::endl;
    std::cout << superc.numberOfSuperbubbles(dfs) << std::endl;
    superc.insert(dfs,g,1384,6936);
    std::cout << "The superbubbles after insertion of an edge are" << std::endl;
    superc.printSuperbubbles(dfs,g);
    std::cout << "The number of superbubbles are" << std::endl;
    std::cout << superc.numberOfSuperbubbles(dfs) << std::endl;
    superc.insert(dfs,g,85,8334);
    std::cout << "The superbubbles after insertion of an edge are" << std::endl;
    superc.printSuperbubbles(dfs,g);
    std::cout << "The number of superbubbles are" << std::endl;
    std::cout << superc.numberOfSuperbubbles(dfs) << std::endl;
    return 0;

}
