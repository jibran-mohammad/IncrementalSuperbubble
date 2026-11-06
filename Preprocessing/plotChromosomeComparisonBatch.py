#!/usr/bin/env python3
"""
plotChromosomeComparisonBatch.py - kept for the old name only.

The comparison plot for both modes is Preprocessing/plotChromosomeComparison.py;
this script calls it in perBatch mode with the same arguments, e.g.

    python3 Preprocessing/plotChromosomeComparisonBatch.py -d ecoli 11 22 X
    = python3 Preprocessing/plotChromosomeComparison.py -m perBatch -d ecoli 11 22 X
"""
import pathlib
import subprocess
import sys

script = pathlib.Path(__file__).resolve().parent / "plotChromosomeComparison.py"
args = sys.argv[1:]
if "-m" not in args and "--mode" not in args:
    args = ["-m", "perBatch"] + args
sys.exit(subprocess.call([sys.executable, str(script)] + args))
