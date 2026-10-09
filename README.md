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
* **g++ / Clang:** A C++ compiler that supports **C++11** or newer (the Makefile compiles with `-std=c++11`).

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
| **minigraph**| The executable must be places directly inside **`Pangenome_Construction/`**. |

## 🚀 Pipeline Execution Guide

The workflow is divided into four phases. You must complete Phase 0 before running any other steps.

### Phase 0: External Tool Setup (One-Time)
This step ensures that the required external bioinformatics tools (vg, mash, minichain, minigraph) are present in the correct locations for the preprocessing scripts to function.

Create the required directory structure:

```Bash

mkdir -p Pangenome_Construction/vgTool
```
**Install/Place Executables**: Install the external tools and ensure their executables are available at the required paths:

**Pangenome_Construction/vgTool/vg**

**Pangenome_Construction/Mash**

**Pangenome_Construction/minichain**

**Pangenome_Construction/minigraph**

### Phase 1: Preprocessing and Graph Generation (graphRunner.py)
This script runs the entire preprocessing pipeline to create the graph file (numeric_edgelist.txt) which is required for Phase 2.

**Script location: scripts/graphRunner.py**

> **Note:** dataset names are case-sensitive and differ between scripts: `graphRunner.py` takes `ecoli` / `Human` / `Virus`, `superbubbleRunner.py` takes `Ecoli` / `Human`, and `parse_to_csv.py` takes `ecoli` / `human` / `all`.

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

**Available Algorithms: `static`, `incRecompute`, `incBST`, `incMinIBST`, `incBSTBatch`.** Each ID is also the source file (`algorithms/<ID>.cpp`), make target, binary (`test/<ID>`) and timing file (`<ID>_time.txt`); see [ALGORITHMS.md](ALGORITHMS.md) for what each one does.

`incBSTBatch` only runs in `perBatch` mode; it is skipped automatically in `perEdge` mode. Each algorithm is run 100 times (`NUM_RUNS` in `superbubbleRunner.py`) and the timings are averaged, so large datasets take a long time.

Examples:
```Bash

# Run all algorithms on E. coli in batch mode
python3 scripts/superbubbleRunner.py --dataset Ecoli -m perBatch

# Run all algorithms on E. coli in Edge mode
python3 scripts/superbubbleRunner.py --dataset Ecoli -m perEdge

# Run specific algorithms on chromosome_1 in batch mode
python3 scripts/superbubbleRunner.py --dataset Human --algorithms incBST --chromosomes chromosome_1 -m perBatch

# Run specific algorithms on chromosome_1 in edge mode
python3 scripts/superbubbleRunner.py --dataset Human --algorithms incBST --chromosomes chromosome_1 -m perEdge
```
### Phase 3: Data Analysis and Plotting
This final phase processes the raw time outputs into readable CSVs and generates performance charts.

#### Step 3a: Generate CSV (parse_to_csv.py)
This script merges the raw time output files from the three algorithms into a single haplotype_build_comparison.csv within the respective timeFile/ directory.

```Bash

# Example: Process E. coli output in Edge mode
python3 Preprocessing/parse_to_csv.py --dataset ecoli -m peredge

# Example: Process E. coli output in Batch mode
python3 Preprocessing/parse_to_csv.py --dataset ecoli -m perbatch

# Example: Process all Human chromosome outputs in Edge mode
python3 Preprocessing/parse_to_csv.py --dataset human -m peredge

# Example: Process all Human chromosome outputs in batch mode
python3 Preprocessing/parse_to_csv.py --dataset human -m perbatch 
```
#### Step 3b: Generate Plots (plotcsv.py)
This script reads the CSV generated in 3a and produces performance charts (PNG files). You must specify the path to the CSV and the desired output directory for the plots.

```Bash

# Example: Generate plots for the E. coli data and save them in a new 'plots/ecoli_charts/' directory in edge mode with Ecoli data set in log scale
python3 Preprocessing/plotcsv.py --csv output/Ecoli/timeFile/perEdge/haplotype_build_comparison.csv --out plots/ecoli_charts/ --mode peredge --dataset Ecoli --scale log

# Example: Generate plots for the E. coli data and save them in a new 'plots/ecoli_charts/' directory in edge mode with Ecoli data set in linear scale
python3 Preprocessing/plotcsv.py --csv output/Ecoli/timeFile/perEdge/haplotype_build_comparison.csv --out plots/ecoli_charts/ --mode peredge --dataset Ecoli --scale linear


# Example: Generate plots for chromosome_11 data in log scale in batch mode.
python3 Preprocessing/plotcsv.py --csv output/Human/chromosome_11/timeFile/perbatch/haplotype_build_comparison.csv --out plots/chr11_charts/ --mode perbatch --dataset GRCh38_chr11 --scale log
```

## 🧪 Test Run Instructions (Quick Start Example)

To verify the installation and test the pipeline end-to-end, use five sample E. coli genomes. Genome files are not tracked in git, so download them from NCBI RefSeq and place the `.fna` files in **input/Ecoli/coli3682_dataset/**:

* `GCF_000005845.2_ASM584v2` (reference)
* `GCF_000007445.1_ASM744v1`
* `GCF_000008865.2_ASM886v2`
* `GCF_000009565.1_ASM956v1`
* `GCF_000010245.2_ASM1024v1`

Follow the steps below to run a complete test execution.

**Step 1** — Ensure Virtual Environment & Dependencies
```Bash
source venv/bin/activate
```
**Step 2** — Preprocessing & Graph Construction

This step processes the five genomes and generates the numeric_edgelist.txt file.
 ```Bash
python3 scripts/graphRunner.py --dataset ecoli
```
**Step 3** — Run Superbubble Algorithms

This will compile the C++ code and run all available algorithms on the generated graph.
```Bash
python3 scripts/superbubbleRunner.py --dataset Ecoli -m perEdge
```
**Step 4** — Generate CSV Output
```Bash
python3 Preprocessing/parse_to_csv.py --dataset ecoli -m peredge
```
**Step 5** — Produce Performance Plots
```Bash
python3 Preprocessing/plotcsv.py --csv output/Ecoli/timeFile/perEdge/haplotype_build_comparison.csv --out output/Ecoli/timeFile --mode peredge --dataset Ecoli --scale log
```
## 📦 Expected Output

After completion, you should find the following:

Location	Contents  
-**input/Ecoli/Graphs/	numeric_edgelist.txt**.  
-**output/Ecoli/timeFile	Algorithm results and timing logs**.  
-**output/Ecoli/timeFile/	Raw timing files + .csv**.  
-**output/Ecoli/timeFile/	Performance comparison graphs (PNG)**.  
 
If everything runs successfully, you have verified the pipeline setup and the algorithms.
