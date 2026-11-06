#include "../lib/graph/contextGraph.h"
#include "../lib/DFS/contextdfs.h"
#include "../lib/findSuperbubble/superbubble.h"
#include "../lib/DFS/simpleDFS.h"
#include "../lib/graph/adjacencyListGraph.h"
#include "../lib/DFS/incrementalDFSAVLTree.h"
#include "../lib/DFS/incrementalDFSNormalBST.h"
#include "../lib/DFS/incrementalDFSAVLTreeH4.h"
#include "../lib/DFS/incrementalDFSSMM.h"
//#include "../lib/DFS/incrementalDFS.h"
#include "../lib/DFS/incDFSBatch.h"
#include "../lib/findSuperbubble/superbubbleh1.h"
#include "../lib/findSuperbubble/superbubbleBatchH1.h"
#include <fstream>
#include <iostream>
#include <string>


/*
  Writes the superbubbles of the final graph next to the timing file:
  <dir>/<ID>_time.txt -> <dir>/<ID>_superbubbles.txt (see ALGORITHMS.md).
  Format: "# Final graph | total=<n>", then one "<entrance> <exit>" per
  line in input ids. Nothing is written if the timing file is not named
  *_time.txt (e.g. /dev/null). Called after all updates, outside the timed
  sections, in normal (non-BENCHMARK) builds only.
*/
template <class SuperbubbleClass>   // SuperbubbleCalculation, SuperbubbleCalculationH1 or SuperbubbleCalculationBatchH1
inline void writeFinalSuperbubbles(const std::string &timePath,
                                   SuperbubbleClass &superbubble,
                                   GraphContext &g) {
    const std::string suffix = "_time.txt";
    if (timePath.size() < suffix.size() ||
        timePath.compare(timePath.size() - suffix.size(), suffix.size(), suffix) != 0) {
        return;
    }
    std::string listPath = timePath.substr(0, timePath.size() - suffix.size()) + "_superbubbles.txt";
    std::ofstream out(listPath);
    if (!out.is_open()) {
        std::cerr << "Error opening superbubble file " << listPath << std::endl;
        return;
    }
    out << "# Final graph | total=" << superbubble.numberOfSuperbubbles() << "\n";
    superbubble.writeSuperbubbles(&g, out);
}


class AlgorithmForSuperbubbles {
public:
     void inputAndCompute();
     long getMemoryUsageInKB();
};

