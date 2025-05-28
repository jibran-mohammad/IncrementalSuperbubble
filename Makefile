# Compiler
CC = g++

# Compiler flags
CFLAGS = -std=c++11 -Wall -g -I../graph/  -I../tree/

# Directories
LIB_DIR = lib
TEST_DIR = test
GRAPH_DIR = $(LIB_DIR)/graph
TREE_DIR = $(LIB_DIR)/tree
DFS_DIR = $(LIB_DIR)/DFS
PREPROCESS_DIR = $(LIB_DIR)/PreprocessingGraphs
RANDOMGRAPH_DIR = $(LIB_DIR)/RandomGraphGeneration
SUPERBUBBLE_DIR = $(LIB_DIR)/findSuperbubble
CYCLEDETECTION_DIR = $(LIB_DIR)/CycleDetection
DAGEXTRACTION_DIR = $(LIB_DIR)/DAGExtraction
COMPARISON_DIR = $(LIB_DIR)/Comparisonalgo
PROBABILITY_DIR = $(LIB_DIR)/ProbabilityEstimation
ALGORITHM_DIR = algorithms
LINKER_FLAGS = -lstdc++

# Source files
GRAPH_SRCS = $(GRAPH_DIR)/adjacencyListGraph.cpp 
CONTEXTGRAPH_SRCS = $(GRAPH_DIR)/contextGraph.cpp
BASEDFS_SRCS = $(DFS_DIR)/baseDFS.cpp
SIMPLEDFS_SRCS = $(DFS_DIR)/simpleDFS.cpp 
DFSCONTEXT_SRCS = $(DFS_DIR)/contextdfs.cpp
INCDFS_SRCS = $(DFS_DIR)/incrementalDFS.cpp
NEWINC_SRCS = $(DFS_DIR)/incrementalDFSAVLTree.cpp
NEWINCDFSH1_SRCS = $(DFS_DIR)/incrementalDFSH1.cpp
SUPERBUBBLE_SRCS = $(SUPERBUBBLE_DIR)/superbubble.cpp
AVLTREE_SRCS = $(TREE_DIR)/avlTree.cpp
DIFFERENTAVLTREE_SRCS = $(TREE_DIR)/differentAVLTree.cpp

# object files
GRAPH_OBJS = $(GRAPH_DIR)/adjacencyListGraph.o $(GRAPH_DIR)/contextGraph.o $(TEST_DIR)/testadjacencyGraph.o
BASEDFS_OBJS = $(DFS_DIR)/baseDFS.o
SIMPLEDFS_OBJS = $(DFS_DIR)/simpleDFS.o $(DFS_DIR)/contextdfs.o $(DFS_DIR)/baseDFS.o $(TEST_DIR)/testSimpleDFS.o
INCDFS_OBJS = $(DFS_DIR)/incremetalDFS.o $(DFS_DIR)/contextdfs.o $(DFS_DIR)/baseDFS.o $(TEST_DIR)/testIncrementalDFS.o
NEWINC_OBJS = $(DFS_DIR)/incrementalDFSAVLTree.o $(DFS_DIR)/contextdfs.o $(DFS_DIR)/baseDFS.o $(TEST_DIR)/testnewIncrementalDFS.o 
NEWINCH1_OBJS = $(DFS_DIR)/incrementalDFSH1.o $(DFS_DIR)/contextdfs.o $(DFS_DIR)/baseDFS.o $(TEST_DIR)/testnewIncrementalDFSH1.o
SUPERBUBBLE_OBJS = $(SUPERBUBBLE_DIR)/superbubble.o 
AVLTREE_OBJS = $(TREE_DIR)/avlTree.o $(TEST_DIR)/testAVLtree.o
DIFFERENTAVLTREE_OBJS = $(TREE_DIR)/differentAVLTree.o $(TEST_DIR)/testDifferentAvlTree.o
SIMPLESUPERBUBBLE_OBJS = $(ALGORITHM_DIR)/superbubbledfs.o
INCSUPERBUBBLE_OBJS = $(ALGORITHM_DIR)/superbubbleIncdfs.o
NEWINCREMENTALSUPERBUBBLE_OBJS = $(ALGORITHM_DIR)/mnlognsuperbubble.o
H1SUPERBUBBLE_OBJS = $(ALGORITHM_DIR)/mnlognH1.o
TESTDFSALLFILES_OBJS = $(TEST_DIR)/testSIIAVLTree.o
TESTALLFILESSUPERBUBBLE_OBJS = $(TEST_DIR)/testAllSuperbubble.o
PREPROCESSGRAPH_OBJS = $(PREPROCESS_DIR)/preprocessGraph.o
GENERATERANDOMGRAPHS_OBJS = $(RANDOMGRAPH_DIR)/generateRandomGraphs.o
GENERATEHUNDEREDRANDOMGRAPHS_OBJS= $(RANDOMGRAPH_DIR)/generateHunderedRandomGraphs.o
isDAG_OBJS= $(CYCLEDETECTION_DIR)/cycleDetection.o
BARANDOMGRAPH_OBJS = $(RANDOMGRAPH_DIR)/barabasiAlbertDAG.o
BARANDOMDIRECTEDGRAPH_OBJS = $(RANDOMGRAPH_DIR)/barabasiAlbertoDG.o
EXTRACTDAGS_OBJS = $(DAGEXTRACTION_DIR)/extractDAGS.o
COMPARISON_OBJS = $(COMPARISON_DIR)/comprisonOfAlgorithms.o
SUPERBUBBLEPROBABILITY_OBJS = $(PROBABILITY_DIR)/randomGraphs.o

# Test source files
TEST_GRAPH_SRCS = $(TEST_DIR)/testadjacencyGraph.cpp
TEST_SIMPLEDFS_SRCS = $(TEST_DIR)/testSimpleDFS.cpp
TEST_INCDFS_SRCS = $(TEST_DIR)/testIncrementalDFS.cpp
TEST_NEWINC_SRCS = $(TEST_DIR)/testnewIncrementalDFS.cpp
TEST_NEWINCH1_SRCS = $(TEST_DIR)/testNewIncrementalDFSh1.cpp
TEST_SUPERBUBBLESIMPLEDFS_SRCS= $(TEST_DIR)/testSuperbubble.cpp
TEST_AVLTREE_SRCS= $(TEST_DIR)/testAVLtree.cpp
TEST_DIFFERENTAVLTREE_SRCS = $(TEST_DIR)/testDifferentAvlTree.cpp
SUPERBUBBLE_SRCS = $(SUPERBUBBLE_DIR)/superbubble.cpp
SIMPLESUPERBUBBLE_SRCS = $(ALGORITHM_DIR)/superbubbledfs.cpp
INCSUPERBUBBLE_SRCS = $(ALGORITHM_DIR)/superbubbleIncdfs.cpp
NEWINCREMENTALSUPERBUBBLE_SRCS = $(ALGORITHM_DIR)/mnlognsuperbubble.cpp
H1SUPERBUBBLE_SRCS = $(ALGORITHM_DIR)/mnlognH1.cpp
TESTDFSALLFILES_SRCS = $(TEST_DIR)/testSIIAVLTree.cpp
TESTALLFILESSUPERBUBBLE_SRCS = $(TEST_DIR)/testAllSuperbubble.cpp
PREPROCESSGRAPH_SRCS = $(PREPROCESS_DIR)/preprocessGraph.cpp
GENERATERANDOMGRAPHS_SRCS = $(RANDOMGRAPH_DIR)/generateRandomGraphs.cpp
GENERATEHUNDEREDRANDOMGRAPHS_SRCS = $(RANDOMGRAPH_DIR)/hunderedRandomGraphs.cpp
isDAG_SRCS= $(CYCLEDETECTION_DIR)/cycleDetection.cpp
BARANDOMGRAPH_SRCS = $(RANDOMGRAPH_DIR)/barabasiAlbertDAG.cpp
BARANDOMDIRECTEDGRAPH_SRCS = $(RANDOMGRAPH_DIR)/barabasiAlbertoDG.cpp
EXTRACTDAGS_SRCS = $(DAGEXTRACTION_DIR)/extractDAGS.cpp
COMPARISON_SRCS = $(COMPARISON_DIR)/comprisonOfAlgorithms.cpp
SUPERBUBBLEPROBABILITY_SRCS = $(PROBABILITY_DIR)/randomGraphs.cpp

# Executables
GRAPH_EXEC = $(TEST_DIR)/testadjacencyListGraph
SIMPLEDFS_EXEC = $(TEST_DIR)/testSimpleDFS
INCDFS_EXEC = $(TEST_DIR)/testIncrementalDFS
NEWINC_EXEC = $(TEST_DIR)/testNewIncrementalDFS
NEWINCH1_EXEC = $(TEST_DIR)/newIncDFSH1
AVLTREE_EXEC = $(TEST_DIR)/testavlTree
DIFFERENTAVLTREE_EXEC = $(TEST_DIR)/differentavlTree
SIMPLESUPERBUBBLE_EXEC = $(TEST_DIR)/superbubbledfs
INCSUPERBUBBLE_EXEC = $(TEST_DIR)/superbubbleIncdfs
NEWINCREMENTALSUPERBUBBLE_EXEC = $(TEST_DIR)/mnlognsuperbubble
H1SUPERBUBBLE_EXEC = $(TEST_DIR)/h1Superbubble
TESTDFSALLFILES_EXEC = $(TEST_DIR)/testAllFilesDfs
TESTALLFILESSUPERBUBBLE_EXEC = $(TEST_DIR)/testAllSuperbubbles
PREPROCESSGRAPH_EXEC = $(PREPROCESS_DIR)/preprocessGraph
GENERATERANDOMGRAPHS_EXEC = $(RANDOMGRAPH_DIR)/generateRandomGraphs
GENERATEHUNDEREDRANDOMGRAPHS_EXEC = $(RANDOMGRAPH_DIR)/generateHunderedRandomGraphs
isDAG_EXEC= $(CYCLEDETECTION_DIR)/cycleDetection
BARANDOMGRAPH_EXEC = $(RANDOMGRAPH_DIR)/barabasiAlbertDAG
BARANDOMDIRECTEDGRAPH_EXEC = $(RANDOMGRAPH_DIR)/barabasiAlbertoDG
EXTRACTDAGS_EXEC = $(DAGEXTRACTION_DIR)/extractdags
COMPARISON_EXEC = $(COMPARISON_DIR)/comparisonAlgorithm
SUPERBUBBLEPROBABILITY_EXEC = $(PROBABILITY_DIR)/randomGraphs

# Default target
# if you run the MAKE command using TARGET=tree, then only tree and its dependencies will
# be handled, if using TARGET=graph, then only graph and its dependencies will be handled 
ifeq ($(TARGET),graph)
	TARGET_EXEC = $(GRAPH_EXEC)
	TARGET_DEPS = $(GRAPH_OBJS) 
else ifeq ($(TARGET),avlTree)
	TARGET_EXEC = $(AVLTREE_EXEC)
	TARGET_DEPS = $(AVLTREE_OBJS)
else ifeq ($(TARGET),differentavlTree)
	TARGET_EXEC = $(DIFFERENTAVLTREE_EXEC)
	TARGET_DEPS = $(DIFFERENTAVLTREE_OBJS)			
else ifeq ($(TARGET),simpledfs)
	TARGET_EXEC = $(SIMPLEDFS_EXEC)
	TARGET_DEPS = $(SIMPLEDFS_OBJS) $(GRAPH_DIR)/adjacencyListGraph.o $(GRAPH_DIR)/contextGraph.o
else ifeq ($(TARGET),incdfs)
	TARGET_EXEC = $(INCDFS_EXEC)
	TARGET_DEPS = $(INCDFS_OBJS) $(GRAPH_DIR)/contextGraph.o $(GRAPH_DIR)/adjacencyListGraph.o $(TREE_DIR)/differentAVLTree.o 
else ifeq ($(TARGET),newincdfs)
	TARGET_EXEC = $(NEWINC_EXEC)
	TARGET_DEPS = $(NEWINC_OBJS) $(GRAPH_DIR)/contextGraph.o $(GRAPH_DIR)/adjacencyListGraph.o $(TREE_DIR)/differentAVLTree.o 
else ifeq ($(TARGET),simpleSuperbubble)
	TARGET_EXEC = $(SIMPLESUPERBUBBLE_EXEC)
	TARGET_DEPS = $(SIMPLESUPERBUBBLE_OBJS) $(SUPERBUBBLE_DIR)/superbubble.o  $(DFS_DIR)/simpleDFS.o $(DFS_DIR)/contextdfs.o $(DFS_DIR)/baseDFS.o $(GRAPH_DIR)/adjacencyListGraph.o	$(GRAPH_DIR)/contextGraph.o
else ifeq ($(TARGET),incsuperbubble)
	TARGET_EXEC = $(INCSUPERBUBBLE_EXEC)
	TARGET_DEPS = $(INCSUPERBUBBLE_OBJS) $(SUPERBUBBLE_DIR)/superbubble.o $(DFS_DIR)/incremetalDFS.o $(DFS_DIR)/contextdfs.o $(DFS_DIR)/baseDFS.o  $(GRAPH_DIR)/adjacencyListGraph.o $(GRAPH_DIR)/contextGraph.o $(TREE_DIR)/differentAVLTree.o
else ifeq ($(TARGET),newincrementalSuperbubble)
	TARGET_EXEC = $(NEWINCREMENTALSUPERBUBBLE_EXEC)
	TARGET_DEPS = $(NEWINCREMENTALSUPERBUBBLE_OBJS) $(SUPERBUBBLE_DIR)/superbubble.o $(DFS_DIR)/incrementalDFSAVLTree.o $(DFS_DIR)/contextdfs.o $(DFS_DIR)/baseDFS.o  $(GRAPH_DIR)/adjacencyListGraph.o $(GRAPH_DIR)/contextGraph.o $(TREE_DIR)/differentAVLTree.o
else ifeq ($(TARGET),H1Superbubble)
	TARGET_EXEC = $(H1SUPERBUBBLE_EXEC)_
	TARGET_DEPS = $(H1SUPERBUBBLE_OBJS) $(SUPERBUBBLE_DIR)/superbubble.o $(DFS_DIR)/incrementalDFSH1.o $(DFS_DIR)/contextdfs.o $(DFS_DIR)/baseDFS.o  $(GRAPH_DIR)/adjacencyListGraph.o $(GRAPH_DIR)/contextGraph.o $(TREE_DIR)/differentAVLTree.o
else ifeq ($(TARGET),testDFSAllFiles)
	TARGET_EXEC = $(TESTDFSALLFILES_EXEC)
	TARGET_DEPS = $(TESTDFSALLFILES_OBJS)
else ifeq ($(TARGET),testAllFilesSuperbubble)
	TARGET_EXEC = $(TESTALLFILESSUPERBUBBLE_EXEC)
	TARGET_DEPS = $(TESTALLFILESSUPERBUBBLE_OBJS)
else ifeq ($(TARGET),preprocessGraph)
	TARGET_EXEC = $(PREPROCESSGRAPH_EXEC)
	TARGET_DEPS = $(PREPROCESSGRAPH_OBJS)
else ifeq ($(TARGET),generateRandomGraphs)
	TARGET_EXEC = $(GENERATERANDOMGRAPHS_EXEC)
	TARGET_DEPS = $(GENERATERANDOMGRAPHS_OBJS)	
else ifeq ($(TARGET),generateHunderedRandomGraphs)
	TARGET_EXEC = $(GENERATEHUNDEREDRANDOMGRAPHS_EXEC)
	TARGET_DEPS = $(GENERATEHUNDEREDRANDOMGRAPHS_OBJS)
else ifeq ($(TARGET),isDAG)
	TARGET_EXEC = $(isDAG_EXEC)
	TARGET_DEPS = $(isDAG_OBJS)	$(GRAPH_DIR)/adjacencyListGraph.o $(GRAPH_DIR)/contextGraph.o
else ifeq ($(TARGET),BARandomGraph)
	TARGET_EXEC = $(BARANDOMGRAPH_EXEC)
	TARGET_DEPS = $(BARANDOMGRAPH_OBJS) $(SUPERBUBBLE_DIR)/superbubble.o  $(DFS_DIR)/simpleDFS.o $(DFS_DIR)/contextdfs.o $(DFS_DIR)/baseDFS.o $(GRAPH_DIR)/adjacencyListGraph.o	$(GRAPH_DIR)/contextGraph.o 
else ifeq ($(TARGET),BARandomDirectedGraph)
	TARGET_EXEC = $(BARANDOMDIRECTEDGRAPH_EXEC)
	TARGET_DEPS = $(BARANDOMDIRECTEDGRAPH_OBJS) 
else ifeq ($(TARGET),extractDAGS)
	TARGET_EXEC = $(EXTRACTDAGS_EXEC)
	TARGET_DEPS = $(EXTRACTDAGS_OBJS) $(GRAPH_DIR)/adjacencyListGraph.o $(GRAPH_DIR)/contextGraph.o	
else ifeq ($(TARGET),comparison)
	TARGET_EXEC = $(COMPARISON_EXEC)
	TARGET_DEPS = $(COMPARISON_OBJS)
else ifeq ($(TARGET),superbubbleProbability)
	TARGET_EXEC = $(SUPERBUBBLEPROBABILITY_EXEC) 
	TARGET_DEPS = $(SUPERBUBBLEPROBABILITY_OBJS) $(SUPERBUBBLE_DIR)/superbubble.o  $(DFS_DIR)/simpleDFS.o $(DFS_DIR)/contextdfs.o $(DFS_DIR)/baseDFS.o $(GRAPH_DIR)/adjacencyListGraph.o	$(GRAPH_DIR)/contextGraph.o
else ifeq ($(TARGET),incdfsH1)
	TARGET_EXEC = $(NEWINCH1_EXEC)
	TARGET_DEPS = $(NEWINCH1_OBJS) $(GRAPH_DIR)/contextGraph.o $(GRAPH_DIR)/adjacencyListGraph.o $(TREE_DIR)/differentAVLTree.o 
endif


# Rules to build target based on argument

# Rule to link executables
$(TARGET_EXEC): $(TARGET_DEPS)
	@echo "Linking $@..."
	@$(CC) $(CFLAGS) $^ -o $@

$(ALGORITHM_DIR)/superbubbledfs.o: $(SIMPLESUPERBUBBLE_SRCS)
	$(CC) $(CFLAGS) -c $< -o $@

$(ALGORITHM_DIR)/superbubbleIncdfs.o: $(INCSUPERBUBBLE_SRCS)
	$(CC) $(CFLAGS) -c $< -o $@	

$(ALGORITHM_DIR)/mnlognsuperbubble.o: $(NEWINCREMENTALSUPERBUBBLE_SRCS)
	$(CC) $(CFLAGS) -c $< -o $@	

$(ALGORITHM_DIR)/mnlognH1.o: $(H1SUPERBUBBLE_SRCS)
	$(CC) $(CFLAGS) -c $< -o $@		

$(SUPERBUBBLE_DIR)/Superbubble.o: $(SUPERBUBBLE_SRCS)
	$(CC) $(CFLAGS) -c $< -o $@
	
$(TESTDFSALLFILES_OBJS): $(TESTDFSALLFILES_SRCS)
	$(CC) $(CFLAGS) -c $< -o $@	

$(TESTALLFILESSUPERBUBBLE_OBJS): $(TESTALLFILESSUPERBUBBLE_SRCS)
	$(CC) $(CFLAGS) -c $< -o $@	

$(GENERATERANDOMGRAPHS_OBJS): $(GENERATERANDOMGRAPHS_SRCS)
	$(CC) $(CFLAGS) -c $< -o $@	

$(GENERATEHUNDEREDRANDOMGRAPHS_OBJS): $(GENERATEHUNDEREDRANDOMGRAPHS_SRCS)
	$(CC) $(CFLAGS) -c $< -o $@		

$(PREPROCESSGRAPH_OBJS): $(PREPROCESSGRAPH_SRCS)
	$(CC) $(CFLAGS) -c $< -o $@	

$(isDAG_OBJS): $(isDAG_SRCS)
	$(CC) $(CFLAGS) -c $< -o $@	

$(BARANDOMGRAPH_OBJS): $(BARANDOMGRAPH_SRCS)
	$(CC) $(CFLAGS) -c $< -o $@	

$(BARANDOMDIRECTEDGRAPH_OBJS): $(BARANDOMDIRECTEDGRAPH_SRCS)
	$(CC) $(CFLAGS) -c $< -o $@	

$(EXTRACTDAGS_OBJS): $(EXTRACTDAGS_SRCS)
	$(CC) $(CFLAGS) -c $< -o $@	

$(COMPARISON_OBJS): $(COMPARISON_SRCS)
	$(CC) $(CFLAGS) -c $< -o $@	

$(SUPERBUBBLEPROBABILITY_OBJS): $(SUPERBUBBLEPROBABILITY_SRCS)
	$(CC) $(CFLAGS) -c $< -o $@					

testIncrementalSuperbubble.o: $(TEST_INCSUPERBUBBLE_SRCS)
	$(CC) $(CFLAGS) -c $< -o $@		

$(SUPERBUBBLE_DIR)/superbubblesimpledfs.o: $(SUPERBUBBLE_SRCS)
	$(CC) $(CFLAGS) -c $< -o $@	

testSuperbubble.o: $(TEST_SUPERBUBBLESIMPLEDFS_SRCS)
	$(CC) $(CFLAGS) -c $< -o $@	
	
$(DFS_DIR)/incDFSNoTree.o: $(INCDFSNOTREE_SRCS)
	$(CC) $(CFLAGS) -c $< -o $@		

$(DFS_DIR)/newIncrementalDFS.o: $(NEWINC_SRCS)
	$(CC) $(CFLAGS) -c $< -o $@	

$(DFS_DIR)/incrementalDFSH1.o: $(NEWINCDFSH1_SRCS)
	$(CC) $(CFLAGS) -c $< -o $@		


testIncrementalDFS.o: $(TEST_INCDFS_SRCS)
	$(CC) $(CFLAGS) -c $< -o $@	

testincDfsNoTree.o: $(TEST_INCDFSNOTREE_SRCS)
	$(CC) $(CFLAGS) -c $< -o $@		

testnewIncrementalDFS.o: $(TEST_NEWINC_SRCS)
	$(CC) $(CFLAGS) -c $< -o $@

	

$(DFS_DIR)/simpleDFS.o: $(SIMPLEDFS_SRCS)
	$(CC) $(CFLAGS) -c $< -o $@

$(DFS_DIR)/incremetalDFS.o: $(INCDFS_SRCS)
	$(CC) $(CFLAGS) -c $< -o $@	

$(DFS_DIR)/baseDFS.o: $(BASEDFS_SRCS)
	$(CC) $(CFLAGS) -c $< -o $@	

$(DFS_DIR)/contextdfs.o: $(DFSCONTEXT_SRCS)
	$(CC) $(CFLAGS) -c $< -o $@	


testSimpleDFS.o: $(TEST_SIMPLEDFS_SRCS)
	$(CC) $(CFLAGS) -c $< -o $@	


adjacencyListTree.o:$(TREE_SRCS)
	$(CC) $(CFLAGS) -I$(TREE_DIR) -c  $< -o $(TREE_DIR)/$@

adjacencyListGraph.o:$(GRAPH_SRCS)
	$(CC) $(CFLAGS) -I$(GRAPH_DIR) -c  $< -o $(GRAPH_DIR)/$@

contextGraph.o: $(CONTEXTGRAPH_SRCS)
	$(CC) $(CFLAGS) -I$(GRAPH_DIR) -c  $< -o $(GRAPH_DIR)/$@	

testadjacencyTree.o:$(TEST_TREE_SRCS)
	$(CC) $(CFLAGS) -I$(TEST_DIR) -c  $< -o $(TEST_DIR)/$@	

testadjacencyGraph.o:$(TEST_GRAPH_SRCS)
	$(CC) $(CFLAGS) -I$(TEST_DIR) -c  $< -o $(TEST_DIR)/$@	

avlTree.o:$(AVLTREE_SRCS)
	$(CC) $(CFLAGS) -I$(TREE_DIR) -c  $< -o $(TREE_DIR)/$@

differentAVLTree.o:$(DIFFERENTAVLTREE_SRCS)
	$(CC) $(CFLAGS) -I$(TREE_DIR) -c  $< -o $(TREE_DIR)/$@
	
testAVLtree.o:$(TEST_AVLTREE_SRCS)
	$(CC) $(CFLAGS) -I$(TREE_DIR) -c  $< -o $(TEST_DIR)/$@

testdifferentAVLtree.o:$(TEST_DIFFERENTAVLTREE_SRCS)
	$(CC) $(CFLAGS) -I$(TREE_DIR) -c  $< -o $(TEST_DIR)/$@	


# Clean rule

clean:
	rm -f $(GRAPH_OBJS) $(GRAPH_EXEC) $(BASEDFS_OBJS) $(SIMPLEDFS_OBJS) $(SIMPLEDFS_EXEC) $(INCDFS_OBJS) $(INCDFS_EXEC) $(NEWINC_OBJS) $(NEWINC_EXEC) $(SUPERBUBBLE_OBJS) $(AVLTREE_OBJS) $(AVLTREE_EXEC) $(DIFFERENTAVLTREE_OBJS) $(DIFFERENTAVLTREE_EXEC) $(SIMPLESUPERBUBBLE_OBJS) $(SIMPLESUPERBUBBLE_EXEC) $(INCSUPERBUBBLE_OBJS) $(INCSUPERBUBBLE_EXEC) $(NEWINCREMENTALSUPERBUBBLE_OBJS) $(NEWINCREMENTALSUPERBUBBLE_EXEC) $(TESTDFSALLFILES_OBJS) $(TESTDFSALLFILES_EXEC) $(TESTALLFILESSUPERBUBBLE_OBJS) $(TESTALLFILESSUPERBUBBLE_EXEC) $(PREPROCESSGRAPH_OBJS) $(PREPROCESSGRAPH_EXEC) $(GENERATERANDOMGRAPHS_OBJS) $(GENERATERANDOMGRAPHS_EXEC) $(GENERATEHUNDEREDRANDOMGRAPHS_OBJS) $(GENERATEHUNDEREDRANDOMGRAPHS_EXEC) $(isDAG_OBJS) $(isDAG_EXEC) $(BARANDOMGRAPH_OBJS) $(BARANDOMGRAPH_EXEC) $(EXTRACTDAGS_OBJS) $(EXTRACTDAGS_EXEC) $(SUPERBUBBLEPROBABILITY_EXEC) $(SUPERBUBBLEPROBABILITY_OBJS) $(NEWINCH1_OBJS) $(NEWINCH1_EXEC)
	unset TARGET
# Phony target to avoid conflicts with files named "clean"
.PHONY: clean
