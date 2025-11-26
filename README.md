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
├── input
├── lib
├── LICENSE
├── Makefile
├── output
├── Pangenome_Construction
├── Preprocessing
├── README.md
├── requirements.txt
├── scripts
├── test
└── venv


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
python3 -m venv venv
source venv/bin/activate

# Install required Python packages
pip install pandas matplotlib numpy
```

### 3. External Executables (CRITICAL)

The pipeline uses several external bioinformatics tools. These executables must be installed and placed in specific locations within the `Pangenome_Construction/` directory before running **Phase 1**.

| Tool | Required Location |
| :--- | :--- |
| **vg** | Must be installed inside the **`Pangenome_Construction/vgTool/`** directory. |
| **mash** | The executable must be placed directly inside **`Pangenome_Construction/`**. |
| **minichain** | The executable must be placed directly inside **`Pangenome_Construction/`**. |

## 🚀 Pipeline Execution Guide

The workflow is divided into four phases. You must complete Phase 0 before running any other steps.

### Phase 0: External Tool Setup (One-Time)
This step ensures that the required external bioinformatics tools (vg, mash, minichain) are present in the correct locations for the preprocessing scripts to function.

Create the required directory structure:

```Bash

mkdir -p Pangenome_Construction/vgTool
```
**Install/Place Executables**: Install the external tools and ensure their executables are available at the required paths:

**Pangenome_Construction/vgTool/vg**

**Pangenome_Construction/Mash**

**Pangenome_Construction/minichain**

### Phase 1: Preprocessing and Graph Generation (graphRunner.py)
This script runs the entire preprocessing pipeline to create the graph file (numeric_edgelist.txt) which is required for Phase 2.

**Script location: scripts/graphRunner.py**

Examples:
```Bash

# Run E. coli Dataset (creates input/Ecoli/Graphs/numeric_edgelist.txt)
python3 scripts/graphRunner.py --dataset ecoli

# Run All Human Chromosomes (iterates through all chromosome_* directories)
python3 scripts/graphRunner.py --dataset Human

# Run Specific Human Chromosomes (e.g., chromosome_1 and chromosome_X)
python3 scripts/graphRunner.py --dataset Human --chromosomes 1 X
```
### Phase 2: Superbubble Algorithm Execution (superbubbleRunner.py)
This step handles C++ compilation via the Makefile and executes the algorithms on the generated graphs, logging raw time output to the output/ directory.

**Script location: scripts/superbubbleRunner.py**

**Available Algorithms: SimpleSuperbubble, IncrementalAlgorithm, AnotherAlgorithm.**

Examples:
```Bash

# Run all algorithms on E. coli
python3 scripts/superbubbleRunner.py --dataset Ecoli

# Run specific algorithms on chromosome_1
python3 scripts/superbubbleRunner.py --dataset Human --algorithms IncrementalAlgorithm --chromosomes chromosome_1
```
### Phase 3: Data Analysis and Plotting
This final phase processes the raw time outputs into readable CSVs and generates performance charts.

#### Step 3a: Generate CSV (parse_to_csv.py)
This script merges the raw time output files from the three algorithms into a single haplotype_build_comparison.csv within the respective timeFile/ directory.

```Bash

# Example: Process E. coli output
python3 Preprocessing/parse_to_csv.py --dataset ecoli

# Example: Process all Human chromosome outputs
python3 Preprocessing/parse_to_csv.py --dataset human 
```
#### Step 3b: Generate Plots (plotcsv.py)
This script reads the CSV generated in 3a and produces performance charts (PNG files). You must specify the path to the CSV and the desired output directory for the plots.

```Bash

# Example: Generate plots for the E. coli data and save them in a new 'plots/ecoli_charts/' directory
python3 plotcsv.py output/Ecoli/timeFile/haplotype_build_comparison.csv plots/ecoli_charts/

# Example: Generate plots for chromosome_1 data
python3 plotcsv.py output/Human/chromosome_1/timeFile/haplotype_build_comparison.csv plots/chr1_charts/
```

## 🧪 Test Run Instructions (Quick Start Example)

To verify the installation and test the pipeline end-to-end, we provide five sample E. coli genomes inside the **input/Ecoli/** directory. Follow the steps below to run a complete test execution.

**Step 1** — Ensure Virtual Environment & Dependencies
```Bash
source venv/bin/activate
```
**Step 2** — Preprocessing & Graph Construction

This step processes the five genomes and generates the numeric_edgelist.txt file.
 ```Bash
python3 scripts/graphRunner.py --dataset Ecoli
```
**Step 3** — Run Superbubble Algorithms

This will compile the C++ code and run all available algorithms on the generated graph.
```Bash
python3 scripts/superbubbleRunner.py --dataset Ecoli
```
**Step 4** — Generate CSV Output
```Bash
python3 Preprocessing/parse_to_csv.py --dataset Ecoli
```
**Step 5** — Produce Performance Plots
python3 Preprocessing/plotcsv.py output/Ecoli/timeFile/haplotype_build_comparison.csv output/Ecoli/timeFile/

## 📦 Expected Output

After completion, you should find the following:

Location	Contents
**input/Ecoli/coli3682_dataset/Graphs/	numeric_edgelist.txt**
**output/Ecoli/timeFile	Algorithm results and timing logs**
**output/Ecoli/timeFile/	Raw timing files + .csv**
**output/Ecoli/timeFile/	Performance comparison graphs (PNG)**

If everything runs successfully, you have verified the pipeline setup and the algorithms.