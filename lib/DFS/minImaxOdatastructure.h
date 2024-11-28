#ifndef MINIMAXODATASTRUCTURE_
#define MINIMAXODATASTRUCTURE_
#include <iostream>
#include <vector>
#include "../graph/contextGraph.h"

class IDatastructure {
public:
	virtual void calculateMinI(GraphContext *g) = 0;
	virtual void calculateMaxO(GraphContext *g) = 0;
	virtual void printMinI() const = 0;
	virtual void printMaxO() const = 0;
	virtual ~IDatastructure() = default;

};
#endif