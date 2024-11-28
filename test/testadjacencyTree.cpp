#include"../lib/tree/tree.h"
#include"../lib/tree/adjacencyListTree.h"


int main()
{
	bool flag;
	AdjacencyListTree adl ;
	Tree *g = &adl;
	g->initialize(15);
	g->addEdge(1,2);
	g->addEdge(2,3);
	g->addEdge(3,4);
	g->addEdge(3,5);
	g->addEdge(3,11);
	g->addEdge(4,8);
	g->addEdge(5,6);
	g->addEdge(5,9);
	g->addEdge(6,7);
	g->addEdge(6,10);
	g->addEdge(8,13);
	g->addEdge(11,12);
	g->addEdge(13,14);
	g->addEdge(13,15);
	g->printTree();
	flag = g->removeEdge(3, 11);
	if(flag == true) {
		std::cout << "the tree after deleting element 11" << std::endl;
		g->printTree();
	} else {
		std::cout << "element not found";
	}
	return 0;


}