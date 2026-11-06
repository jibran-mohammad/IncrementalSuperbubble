#ifndef DFS_SIMPLE
#define DFS_SIMPLE
#include "baseDFS.h"
#include "minImaxOdatastructure.h"
#include <climits>

class SimpleDFS : public BaseDFS , public IDatastructure {
public:
	SimpleDFS(GraphContext *g): BaseDFS(g) {
    }
	// true (perEdge): every insert clears the arrays for the DFS that follows
	// the edge (BaseDFS::insert). false (perBatch): insert only adds the edge
	// and grows the arrays; the program calls clear() once before the DFS of
	// the batch.
	bool clearOnInsert = true;
		
public:
	void insert(GraphContext *g, int u, int v) override;
	void clear();
	void calculateMinI(GraphContext *g) override;
	void calculateMaxO(GraphContext *g) override;
	void printMinI() const override;
	void printMaxO() const override;
};
#endif