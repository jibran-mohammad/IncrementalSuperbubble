"""
bubblefinderTool.py

BubbleFinder implementation of the superbubble tool interface
(toolInterface.py), used by scripts/superbubbleComparator.py -t bubblefinder.

BubbleFinder (Sena, Politov et al., arXiv:2511.21919) computes superbubbles
in linear time with SPQR trees. It is a separate program (C++/Rust), run
once per update on the same temporary GFA that BubbleGun gets:

    BubbleFinder superbubbles --directed --gfa-directed -j 1 \
        -g <gfa> -o <out> --report-json <json>

    --directed --gfa-directed   the graph is directed (our links are all + +);
                                checked on the worked example
                                (examples/regionAnalysis, branch region-analysis)
                                and on chromosome 21: the same superbubbles as
                                the reference / Inc-MinI-BST / BubbleGun
    -j 1                        one thread, like our algorithms and BubbleGun

=====================================================================
WHAT IS MEASURED
=====================================================================

time    BubbleFinder's own timers from --report-json: reading the GFA
        (io/read_graph) + every phase of finding the superbubbles (sb/...).
        The same work as BubbleGun's "time" (reading the GFA + find_bubbles):
        program start-up (about 10 ms) and writing the output file are not
        included.
phases  read, find (as above) and write (writing the output, not in "time")
memory  peak resident memory of the BubbleFinder process (VmHWM, the same
        measure as PeakMemory(KB) of the C++ algorithms), from the report of
        a separate, untimed run; baseline = the peak of BubbleFinder on a
        graph with one edge (program and libraries without a graph)

Bubbles: every superbubble BubbleFinder reports, as ("superbubble", entrance,
exit) with GFA segment names; it writes the entrance first and adds an
orientation sign ("s12+"), which is removed. BubbleFinder reports all
superbubbles, also those with a single edge (unlike BubbleGun); the agreement
table leaves single-edge superbubbles out on both sides anyway.

=====================================================================
INSTALLATION
=====================================================================

Linux binary from https://github.com/algbio/BubbleFinder/releases (v1.0.5):

    mkdir -p ~/tools/bubblefinder-1.0.5 && cd ~/tools/bubblefinder-1.0.5
    wget https://github.com/algbio/BubbleFinder/releases/download/v1.0.5/BubbleFinder-linux-x86_64.tar.gz
    tar xzf BubbleFinder-linux-x86_64.tar.gz
    ./BubbleFinder --version

Default path: ~/tools/bubblefinder-1.0.5/BubbleFinder; another one with the
environment variable BUBBLEFINDER=/path/to/BubbleFinder.
"""

import json
import os
import shutil
import subprocess
import tempfile

from .toolInterface import SuperbubbleTool


DEFAULT_BINARY = os.path.expanduser("~/tools/bubblefinder-1.0.5/BubbleFinder")


class BubbleFinderTool(SuperbubbleTool):

    name = "bubblefinder"

    def __init__(self):

        self.binary = os.environ.get("BUBBLEFINDER", DEFAULT_BINARY)

        if not (os.path.isfile(self.binary) and os.access(self.binary, os.X_OK)):
            raise RuntimeError(
                f"BubbleFinder not found at {self.binary}; install it (see "
                f"Preprocessing/superbubbleTools/bubblefinderTool.py) or set BUBBLEFINDER")

        # Private folder for the output and report of each run
        self.tmp = tempfile.mkdtemp(prefix="bubblefinder_")
        self.out = os.path.join(self.tmp, "superbubbles.txt")
        self.report = os.path.join(self.tmp, "report.json")
        self.baseline = None


    def _call(self, gfa_file):
        """One BubbleFinder run; returns the marks of its JSON report."""

        for f in (self.out, self.report):
            if os.path.exists(f):
                os.remove(f)

        result = subprocess.run(
            [self.binary, "superbubbles", "--directed", "--gfa-directed", "-j", "1",
             "-g", gfa_file, "-o", self.out, "--report-json", self.report],
            capture_output=True, text=True)

        if result.returncode != 0:
            tail = "\n".join((result.stderr or result.stdout).strip().splitlines()[-5:])
            raise RuntimeError(f"BubbleFinder failed on {gfa_file} (return code {result.returncode}):\n{tail}")

        return json.load(open(self.report))["marks"]


    def run(self, gfa_file, keep_bubbles=False):

        marks = self._call(gfa_file)

        read = sum(m["seconds"] for m in marks if m["label"] == "io/read_graph")
        find = sum(m["seconds"] for m in marks if m["label"].startswith("sb/"))
        write = sum(m["seconds"] for m in marks if m["label"] == "io/write_output")

        # Output: first line the number of superbubbles, then "<entrance> <exit>"
        lines = open(self.out).read().split("\n")
        n = int(lines[0])
        pairs = [line.split() for line in lines[1:] if line.strip()]

        if len(pairs) != n:
            raise RuntimeError(f"BubbleFinder output of {gfa_file}: {n} announced, {len(pairs)} listed")

        bubbles = []
        if keep_bubbles:
            bubbles = [("superbubble", a.rstrip("+-"), b.rstrip("+-")) for a, b in pairs]

        return {
            "time": read + find,
            "phases": {"read": read, "find": find, "write": write},
            "counts": {"superbubble": n},
            "bubbles": bubbles,
        }


    def memory(self, gfa_file):

        # Baseline once: BubbleFinder on a graph with a single edge
        if self.baseline is None:
            tiny = os.path.join(self.tmp, "tiny.gfa")
            with open(tiny, "w") as f:
                f.write("H\tVN:Z:1.0\nS\ta\t*\nS\tb\t*\nL\ta\t+\tb\t+\t0M\n")
            self.baseline = max(m["hwm_max_bytes"] for m in self._call(tiny)) // 1024

        peak = max(m["hwm_max_bytes"] for m in self._call(gfa_file)) // 1024

        return {"peak": peak, "baseline": self.baseline}


    def close(self):

        shutil.rmtree(self.tmp, ignore_errors=True)
