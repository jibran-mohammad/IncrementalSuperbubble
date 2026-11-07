# IncSuperbubble: Incremental Superbubble Maintenance in Pangenome Graphs

This repository contains the C++ algorithms that keep the superbubbles of a
growing pangenome DAG up to date with an incremental DFS, and the Python
scripts that run the experiments of the paper.

The algorithms, with the names used in the paper:

| Name | Per edge | Per batch (haplotype) |
|---|---|---|
| Simple | `static` | `static` |
| ImpDFS | `incRecompute` | `incRecomputeBatch` |
| ImpDS | `incMinIBST` | `incBSTBatch` |
| ImpSB | `incMinIBSTH1` | `incBSTBatchH2` |

Each ID is the source file (`algorithms/<ID>.cpp`), the make target and the
binary (`test/<ID>`).

---

## 🧪 Testing the Artifact

Two scripts reproduce the experiments of the paper. Run them from the
repository root.

### Getting the code

This needs `git` (installation: see Requirements below).

```Bash
git clone https://github.com/jibran-mohammad/IncrementalSuperbubble.git
cd IncrementalSuperbubble
```

Every command below is run from this folder. The test data (`testdata/`)
comes with the clone; nothing else needs to be downloaded for test 1.

### Requirements

* Linux (x86-64), `git`, `make`, `g++` (C++11) and Python 3.8 or newer
  with `venv` and `pip`. Check with `git --version`, `make --version`,
  `g++ --version` and `python3 --version`. To install what is missing (needs
  `sudo`):

  ```Bash
  # Debian / Ubuntu
  sudo apt-get update && sudo apt-get install -y git make g++ python3 python3-venv python3-pip

  # Fedora / RHEL / Rocky
  sudo dnf install -y git make gcc-c++ python3 python3-pip
  ```
* The Python packages, including BubbleGun:

  ```Bash
  python3 -m venv venv && source venv/bin/activate
  pip install -r requirements.txt
  ```
* BubbleFinder 1.0.5. The scripts look for it at
  `~/tools/bubblefinder-1.0.5/BubbleFinder`, or set
  `BUBBLEFINDER=/path/to/BubbleFinder`:

  ```Bash
  mkdir -p ~/tools/bubblefinder-1.0.5 && cd ~/tools/bubblefinder-1.0.5 \
    && wget https://github.com/algbio/BubbleFinder/releases/download/v1.0.5/BubbleFinder-linux-x86_64.tar.gz \
    && tar xzf BubbleFinder-linux-x86_64.tar.gz && cd -
  ```
* For test 2 only: `curl`, `vg`, about 10 GB of free disk space and 8 GB
  of memory.
  * `curl` (downloads the HPRC graph). Check with `curl --version`; if it
    is missing, install it (Debian/Ubuntu shown; it needs `sudo`):

    ```Bash
    sudo apt-get update && sudo apt-get install -y curl
    ```
  * `vg` (converts the HPRC graph to GFA): a single static Linux binary
    from https://github.com/vgteam/vg/releases. No installation or `sudo`
    is needed:

    ```Bash
    mkdir -p ~/tools && cd ~/tools \
      && wget https://github.com/vgteam/vg/releases/download/v1.78.0/vg \
      && chmod +x vg && ./vg version && cd -
    ```

    Then give it to the test with `--vg ~/tools/vg`.

### Test 1: the small human graphs (`scripts/smallGraphsTest.py`)

The five small graphs of the paper (chromosomes 11, 17, 18, 22 and X) are
in `testdata/small/` (gzipped; base graph, haplotypes and the edge list of
the updates). The script unpacks them and then:

* runs **our algorithms** per edge and per batch (the four above, and per
  edge also the data structures NBST, `incNormalBST`, and BST, `incBST`);
* runs **BubbleGun** and **BubbleFinder** per batch;
* prints the **tables** and writes them to `output/artifactTest/small_tables.md`.
  They show the total time per graph with the speedup over Simple, the time
  of every step, minI/maxO for each data structure, the tools against our
  algorithms (time and memory), and whether all our algorithms and the
  tools find the same superbubbles in the final graph. The time of the tools
  is also split into reading the GFA and the search alone, because they
  re-read the whole graph after every update.

```Bash
python3 scripts/smallGraphsTest.py            # all five graphs: about 40 minutes
python3 scripts/smallGraphsTest.py -c 22      # one graph only: a few minutes
```

**Per edge, the tools can take hours.** BubbleGun and BubbleFinder cannot
update a result. After every inserted edge they read and search the whole
graph again, about 60,000 times for the five graphs, so this can take many
hours, depending on the machine. Per batch they run once per haplotype, a
few hundred times. To also run the tools per edge (optional):

```Bash
python3 scripts/smallGraphsTest.py --tools-per-edge
```

Other options: `--runs N` (runs of each of our algorithms, averaged; default
1, the paper used 30) and `--no-tools`. Absolute times depend on the machine;
the speedups much less.

### Test 2: the conversion of an HPRC graph into a DAG (`scripts/dagConversionTest.py`)

This test checks how much the conversion of a Minigraph-Cactus graph into a
DAG (`Preprocessing/mcToDag.py`) changes the superbubbles. It takes one
chromosome of HPRC release 1.1 (44 samples, on GRCh38). The graph is too
large for this repository, so the script downloads it from the HPRC bucket.
The script then:

* converts the graph into a DAG and reports what was removed;
* finds the superbubbles of the original graph with BubbleFinder
  (bidirected), and those of the DAG with our algorithm ImpSB after the last
  haplotype;
* as a check, also runs BubbleFinder on the DAG, which must give the same
  superbubbles as ImpSB;
* writes the tables to `output/artifactTest/dag_conversion_chr<c>.md`.

Choose the chromosome with `-c`:

| `-c` | Download | Superbubbles of the original lost in the DAG |
|---|---:|---:|
| 22 | 0.94 GB | 1.16 % |
| 18 | 1.51 GB | 0.00 % |
| 17 | 1.54 GB | 0.82 % |
| X  | 1.83 GB | 0.19 % |
| 11 | 2.50 GB | 0.20 % |

```Bash
python3 scripts/dagConversionTest.py -c 22 --vg ~/tools/vg
```

The download runs at about 1 MB/s (chr22: about 15 minutes); the rest takes
a few minutes. Every other chromosome (1–22, X, Y) works the same way.

### Data

The small graphs were built with minigraph/minichain from assemblies of the
Human Pangenome Reference Consortium (HPRC, Data Release 2) on the
T2T-CHM13 reference. The HPRC graphs are the consortium's Minigraph-Cactus
graphs (release 1.1 and 2.1, on GRCh38). We thank the HPRC for making these
data publicly available.

---

## 📁 Repository Layout

```
algorithms/      one C++ file per algorithm (main program)
lib/             incremental DFS, data structures, superbubble computation
Makefile         make TARGET=<ID> builds test/<ID>
scripts/         run scripts: the two tests above, superbubbleRunner.py
                 (our algorithms), superbubbleComparator.py (the tools),
                 hprcPrepare.py (HPRC download and conversion)
Preprocessing/   graph preparation (mcToDag.py), parsing, tables, plots
testdata/small/  the five small graphs
```

Every script has a guide at the top of its file, and
`python3 <script> --help` lists its options.

## License

GNU General Public License v3.0, see [LICENSE](LICENSE).
