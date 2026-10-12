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

Every script you run has a guide at the top of its file (HOW TO RUN, options, input and output files), and `python3 <script> --help` lists its options. Algorithm names are explained in [ALGORITHMS.md](ALGORITHMS.md); how the graphs, base graphs and updates are prepared (rules, checks, dataset numbers) is described in [GRAPH_PREPARATION.md](GRAPH_PREPARATION.md); how the algorithm is compared with external tools such as BubbleGun (what is timed, memory, agreement, how to read the speedup, which older results are not comparable) in [COMPARISON_METHOD.md](COMPARISON_METHOD.md).

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

#### Step 3b2: Check the superbubbles after every update (checkSuperbubbles.py)
Compares the superbubble count of one or more algorithms after every update with an independent reference (and the final superbubble lists), e.g.:

```Bash
python3 Preprocessing/checkSuperbubbles.py --edgelist input/Ecoli/Graphs/numeric_edgelist.txt -m perEdge \
    --time output/Ecoli/timeFile/perEdge/static_time.txt output/Ecoli/timeFile/perEdge/incMinIBST_time.txt
```

#### Step 3c: Compare the best algorithm with an external tool (toolComparisonToLatex.py)
Writes LaTeX tables comparing the best algorithm (`BEST_ALGORITHM` in `scripts/pipelineRunner.py`) with a superbubble tool (default `bubblegun`; any tool run through `scripts/superbubbleComparator.py -t <tool>` works the same way). It needs, per dataset and mode:

* the tool results: `python3 scripts/superbubbleComparator.py -d Ecoli -m perEdge`
* the algorithm time: `python3 scripts/pipelineRunner.py -d Ecoli -m perEdge`
* the algorithm memory and superbubble list: `python3 scripts/superbubbleRunner.py -d Ecoli -m perEdge`

```Bash
python3 Preprocessing/toolComparisonToLatex.py -d Ecoli
python3 Preprocessing/toolComparisonToLatex.py -d Ecoli Human -c 11 13 X -t bubblegun
```
Output in `output/toolComparison/<tool>/`: `performance_<mode>.tex` (time, speedup, peak memory), `agreement.tex` (superbubbles of the final graph found by both / only one) and `agreement_details_<dataset>.txt` (the differing superbubbles).

## 🔁 Rerunning All Experiments (e.g. on a Server)

Use this after a change to the preprocessing or to the algorithms. The
pangenome graphs (`<ref>.gfa` with walks) do **not** have to be generated
again: everything below starts from the existing graphs. Do not run
`scripts/graphRunner.py` for this (it rebuilds the graphs, which takes
months for the human chromosomes).

**0. Code and data**

```Bash
git pull                      # branch repo-new-track
source venv/bin/activate
pip install -r requirements.txt
```

Every dataset needs its graph with walks:
`input/Ecoli/Graphs/GCF_000005845.2_ASM584v2_genomic.gfa` and
`input/Human/chromosome_<X>/Graphs/chromosome_<X>.gfa`. The next step
rewrites the other files in these `Graphs/` folders; keep a copy if you
want the old ones.

**1. Preprocessing + end-to-end time of the best algorithm** (`pipelineRunner.py`)

Regenerates `reduced_graph.gfa`, `unique_elements.txt`,
`ordered_haplotypes.txt`, `numeric_edgelist.txt` and
`haplotype_removal_report.txt` (rules and checks:
[GRAPH_PREPARATION.md](GRAPH_PREPARATION.md)). One command does this for
E. coli and every human chromosome that has a graph:

```Bash
python3 scripts/pipelineRunner.py -d all --preprocess-only
```

It takes seconds per chromosome. A chromosome that fails a check is
skipped and the others continue; the SUMMARY at the end lists every
dataset with its status and graph sizes, and the command exits with
code 1 if any failed. Then time the best algorithm (10 runs) for both
modes:

```Bash
python3 scripts/pipelineRunner.py -d all -m perEdge 2>&1 | tee pipeline_perEdge.log
python3 scripts/pipelineRunner.py -d all -m perBatch 2>&1 | tee pipeline_perBatch.log
```

(These also rerun the preprocessing, with the same result.)

Chromosomes without a graph are skipped. Each chromosome's
counts are in `input/Human/chromosome_<X>/Graphs/haplotype_removal_report.txt`.

**2. All algorithms** (`superbubbleRunner.py`, 100 runs each): timing
files, peak memory and final superbubble lists for the plots and tables.
Only one `superbubbleRunner.py` at a time per copy of the repository.

```Bash
python3 scripts/superbubbleRunner.py -d Ecoli -m perEdge
python3 scripts/superbubbleRunner.py -d Ecoli -m perBatch -a static incRecompute incBST incMinIBST
python3 scripts/superbubbleRunner.py -d Human -m perEdge
python3 scripts/superbubbleRunner.py -d Human -m perBatch -a static incRecompute incBST incMinIBST
```

(Leave out `-a ...` in perBatch once `incBSTBatch` compiles.)

**3. BubbleGun** (`superbubbleComparator.py`):

```Bash
python3 scripts/superbubbleComparator.py -d all -m perEdge  2>&1 | tee bubblegun_perEdge.log
python3 scripts/superbubbleComparator.py -d all -m perBatch 2>&1 | tee bubblegun_perBatch.log
```

`-d all` runs E. coli and every preprocessed human chromosome one after
the other and ends with a SUMMARY (`-d Human -c 14 X` for single
chromosomes).

**Rough times** (measured on chromosome 21: 64 genomes, 6,863 edge
updates; larger chromosomes take longer): step 1 about 1 minute per
chromosome; step 2 about 50 minutes perEdge (half of it for `static`) and
7 minutes perBatch; step 3 about 30 seconds perBatch and, by estimate,
30-60 minutes perEdge.

**4. Copy back** (to the same paths in this repository):

```
output/Ecoli/{timeFile,incsuperbubble,bubblegun}/
output/Human/chromosome_*/{timeFile,incsuperbubble,bubblegun}/
input/Ecoli/Graphs/{numeric_edgelist.txt,haplotype_removal_report.txt}
input/Human/chromosome_*/Graphs/{numeric_edgelist.txt,haplotype_removal_report.txt}
```

(`numeric_edgelist.txt` is needed by `toolComparisonToLatex.py`; the
reports give the graph sizes for the paper.)

**5. Tables and plots** (locally): Phase 3 above
(`parse_to_csv.py -d all`, `plotcsv.py`, `csvToLatex.py`) and Step 3c
(`toolComparisonToLatex.py -d Ecoli Human`).

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
