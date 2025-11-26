# IncSuperbubble: Pangenome Superbubble Detection Pipeline

This repository contains the source code for the **IncSuperbubble** C++ algorithms and the Python scripts used to orchestrate the entire pangenome processing pipeline.

The pipeline is generally run in two phases:
1.  **Preprocessing** (handled by `graphRunner.py`): Converts raw genome data into processed, normalized, and simplified edge list files (`numeric_edgelist.txt`).
2.  **Superbubble Detection** (handled by `superbubbleRunner.py`): Compiles and executes the C++ superbubble algorithms on the generated edge lists, measuring their performance.

---

## 🏗️ Directory Structure

The project uses the following directory structure, rooted at `IncSuperbubble/`:
```bash
.
├── algorithms
│   ├── mnlognsuperbubble.cpp
│   ├── mnlognsuperbubble.o
│   ├── superbubbleAlgorithms.h
│   ├── superbubbledfs.cpp
│   ├── superbubbledfs.o
│   ├── superbubbleIncdfs.cpp
│   └── superbubbleIncdfs.o
├── input
│   ├── Ecoli
│   ├── HaplotypeTestingData
│   └── Human
├── lib
│   ├── Camparisonalgo
│   ├── CycleDetection
│   ├── DAGExtraction
│   ├── DFS
│   ├── findSuperbubble
│   ├── graph
│   ├── plot
│   ├── PreprocessingGraphs
│   ├── ProbabilityEstimation
│   ├── RandomGraphGeneration
│   └── tree
├── LICENSE
├── Makefile
├── output
│   ├── Ecoli
│   ├── graphs_incremental
│   ├── Human
│   ├── performance_plot_GRC.png
│   ├── testDAG.gfa
│   └── Yeast
├── Pangenome_Construction
│   ├── agc
│   ├── bifrost
│   ├── home
│   ├── Mash
│   ├── minichain
│   ├── minigraph
│   ├── odgi
│   └── vgTool
├── Preprocessing
│   ├── csvToLatex.py
│   ├── DAGCheck.py
│   ├── differencesInGraphs.py
│   ├── discover_hprc_urls.sh
│   ├── download_chromosomes_https.sh
│   ├── duplicateSegments.py
│   ├── edgeListProcessor.py
│   ├── fetch_ecoli_genomes.sh
│   ├── gfaToEdgeListConverter.py
│   ├── graphgeneration.py
│   ├── haplotypeRemoval.py
│   ├── haplotypeSorter.py
│   ├── inconsistencyCheck.py
│   ├── normalizeOrientations.py
│   ├── organize_chromosomes.py
│   ├── parse_to_csv.py
│   ├── plotcsv.py
│   ├── reverseComplement.py
│   ├── split_fasta.py
│   ├── subsetCheck.py
│   ├── test_path.py
│   └── typesOfLinks.py
├── README.md
├── requirements.txt
├── scripts
│   ├── graphRunner.py
│   └── superbubbleRunner.py
├── test
│   ├── changedgraph.txt
│   ├── differentgraph.txt
│   ├── graph1.txt
│   ├── graph2.txt
│   ├── graphProcessed.txt
│   ├── graphsir.txt
│   ├── graph.txt
│   ├── h1Superbubble_
│   ├── mnlognsuperbubble
│   ├── newgraph.txt
│   ├── randomgraph2.txt
│   ├── superbubbledfs
│   ├── superbubbleIncdfs
│   ├── testadjacencyGraph.cpp
│   ├── testadjacencyTree.cpp
│   ├── testAllSuperbubble.cpp
│   ├── testAVLtree.cpp
│   ├── testDifferentAvlTree.cpp
│   ├── testgraph.txt
│   ├── testIncrementalDFS.cpp
│   ├── testnewIncrementalDFS.cpp
│   ├── testnewIncrementalDFSH1.cpp
│   ├── testSIIAVLTree.cpp
│   ├── testSimpleDFS.cpp
│   └── testSuperbubble.cpp
└── venv
    ├── bin
    ├── include
    ├── lib
    ├── lib64 -> lib
    ├── pyvenv.cfg
    └── share

```
---

## 🛠️ Prerequisites

### 1. General System Tools
* **Git:** For cloning the repository.
* **Python 3:** Version 3.6 or newer.
* **make:** The build automation tool.
* **g++ / Clang:** A C++ compiler that supports **C++17** or newer.

### 2. Python Dependencies
The following libraries are required for the pipeline, especially for the analysis and plotting phases.

```bash
# Recommended: Create and activate a virtual environment first
# python3 -m venv venv
# source venv/bin/activate

# Install required Python packages
pip install pandas matplotlib numpy

### 3. External Executables (CRITICAL)

The pipeline uses several external bioinformatics tools. These executables must be installed and placed in specific locations within the `Pangenome_Construction/` directory before running **Phase 1**.

| Tool | Required Location |
| :--- | :--- |
| **vg** | Must be installed inside the **`Pangenome_Construction/vgTool/`** directory. |
| **mash** | The executable must be placed directly inside **`Pangenome_Construction/`**. |
| **minichain** | The executable must be placed directly inside **`Pangenome_Construction/`**. |

🚀 Execution Instructions
The pipeline must be run sequentially in four distinct phases. You must complete Phase 0 before starting any other phase.

Phase 0: External Tool Setup (One-Time)
This step ensures that the required external bioinformatics tools (vg, mash, minichain) are present in the correct locations for the preprocessing scripts to function.

Create the required directory structure:

Bash

mkdir -p Pangenome_Construction/vgTool
Install/Place Executables: Install the external tools and ensure their executables are available at the required paths:

Pangenome_Construction/vgTool/vg

Pangenome_Construction/Mash

Pangenome_Construction/minichain

Phase 1: Preprocessing and Graph Generation (graphRunner.py)
This script runs the entire preprocessing pipeline to create the graph file (numeric_edgelist.txt) which is required for Phase 2.

Script location: scripts/graphRunner.py

Examples:
Bash

# Run E. coli Dataset (creates input/Ecoli/Graphs/numeric_edgelist.txt)
python3 scripts/graphRunner.py --dataset ecoli

# Run All Human Chromosomes (iterates through all chromosome_* directories)
python3 scripts/graphRunner.py --dataset Human

# Run Specific Human Chromosomes (e.g., chromosome_1 and chromosome_X)
python3 scripts/graphRunner.py --dataset Human --chromosomes 1 X
Phase 2: Superbubble Algorithm Execution (superbubbleRunner.py)
This step handles C++ compilation via the Makefile and executes the algorithms on the generated graphs, logging raw time output to the output/ directory.

Script location: scripts/superbubbleRunner.py

Available Algorithms: SimpleSuperbubble, IncrementalAlgorithm, AnotherAlgorithm.

Examples:
Bash

# Run all algorithms on E. coli
python3 scripts/superbubbleRunner.py --dataset Ecoli

# Run specific algorithms on chromosome_1
python3 scripts/superbubbleRunner.py --dataset Human --algorithms IncrementalAlgorithm --chromosomes chromosome_1
Phase 3: Data Analysis and Plotting
This final phase processes the raw time outputs into readable CSVs and generates performance charts.

Step 3a: Generate CSV (parse_to_csv.py)
This script merges the raw time output files from the three algorithms into a single haplotype_build_comparison.csv within the respective timeFile/ directory.

Bash

# Example: Process E. coli output
python3 Preprocessing/parse_to_csv.py --dataset ecoli

# Example: Process all Human chromosome outputs
python3 Preprocessing/parse_to_csv.py --dataset human 
Step 3b: Generate Plots (plotcsv.py)
This script reads the CSV generated in 3a and produces performance charts (PNG files). You must specify the path to the CSV and the desired output directory for the plots.

Bash

# Example: Generate plots for the E. coli data and save them in a new 'plots/ecoli_charts/' directory
python3 plotcsv.py output/Ecoli/timeFile/haplotype_build_comparison.csv plots/ecoli_charts/

# Example: Generate plots for chromosome_1 data
python3 plotcsv.py output/Human/chromosome_1/timeFile/haplotype_build_comparison.csv plots/chr1_charts/

