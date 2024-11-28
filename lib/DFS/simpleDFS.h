#ifndef DFS_SIMPLE
#define DFS_SIMPLE
#include "baseDFS.h"
#include "minImaxOdatastructure.h"
#include <climits>

class SimpleDFS : public BaseDFS , public IDatastructure {
public:
	SimpleDFS(GraphContext *g): BaseDFS(g) {
    }
		
public:
	void calculateMinI(GraphContext *g) override;
	void calculateMaxO(GraphContext *g) override;
	void printMinI() const override;
	void printMaxO() const override;
};
#endif