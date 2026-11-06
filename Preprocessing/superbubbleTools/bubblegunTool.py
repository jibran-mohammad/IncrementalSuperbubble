"""
bubblegunTool.py

BubbleGun implementation of the superbubble tool interface.

BubbleGun is used as a Python library in this process instead of being
launched as a program for every update, so the measured time does not
include starting Python and importing BubbleGun (about 50 ms per run,
roughly half of a whole BubbleGun run on the E. coli graph).

Timed phases of one run (seconds):

    read      parse the whole GFA into BubbleGun's graph
    find      find_bubbles(): detect all bubbles from scratch
    chains    connect bubbles into chains and find nesting
              (extra BubbleGun work, not part of superbubble detection)

"time" = read + find: what a static tool has to redo after every update.

Memory (memory()) is measured in a fresh Python process, because the peak
memory of this process only ever grows and includes the comparator itself:
peak = VmHWM after read + find (the same work as "time"), baseline = VmHWM
after importing BubbleGun, before reading the graph.

BubbleGun reports three bubble types (simple, super, insertion);
all of them are superbubbles, so all three are counted.
"""

import logging
import subprocess
import sys
import time

from .toolInterface import SuperbubbleTool

from BubbleGun.Graph import Graph
from BubbleGun.find_bubbles import find_bubbles
from BubbleGun.connect_bubbles import connect_bubbles
from BubbleGun.find_parents import find_parents


# Run in a fresh interpreter by BubbleGunTool.memory(); prints
# "<baseline KB> <peak KB>" (VmHWM = peak resident memory of this process).
MEMORY_SCRIPT = """
import logging, sys

def vm_hwm():
    for line in open("/proc/self/status"):
        if line.startswith("VmHWM:"):
            return int(line.split()[1])

logging.disable(logging.CRITICAL)
from BubbleGun.Graph import Graph
from BubbleGun.find_bubbles import find_bubbles

baseline = vm_hwm()
graph = Graph(sys.argv[1])
find_bubbles(graph)
print(baseline, vm_hwm())
"""


class WarningCounter(logging.Handler):
    """Counts BubbleGun warnings (e.g. links to undefined segments)."""

    def __init__(self):

        super().__init__(level=logging.WARNING)
        self.count = 0
        self.first = None

    def emit(self, record):

        self.count += 1

        if self.first is None:
            self.first = record.getMessage()


class BubbleGunTool(SuperbubbleTool):

    name = "bubblegun"

    def __init__(self):

        # BubbleGun logs through the root logger; keep its messages out of
        # the terminal and count the warnings instead.
        self.warnings = WarningCounter()

        root = logging.getLogger()
        root.addHandler(self.warnings)
        root.setLevel(logging.WARNING)


    def run(self,gfa_file,keep_bubbles=False):

        self.warnings.count = 0
        self.warnings.first = None

        t0=time.perf_counter()

        graph=Graph(gfa_file)

        t1=time.perf_counter()

        find_bubbles(graph)

        t2=time.perf_counter()

        connect_bubbles(graph)
        find_parents(graph)

        t3=time.perf_counter()

        simple,super_,insertion=graph.bubble_number()

        bubbles=[]

        # Same type names and ends as BubbleGun's --bubble_json output.
        if keep_bubbles:

            for chain in graph.b_chains:
                for bubble in chain.bubbles:

                    if bubble.is_simple():
                        bubble_type="simple"
                    elif bubble.is_super():
                        bubble_type="super"
                    else:
                        bubble_type="insertion"

                    bubbles.append((bubble_type,bubble.source.id,bubble.sink.id))

        if self.warnings.count:
            raise RuntimeError(
                f"BubbleGun reported {self.warnings.count} warning(s) on {gfa_file}, "
                f"first: {self.warnings.first}"
            )

        return {
            "time": (t1-t0)+(t2-t1),
            "phases": {
                "read": t1-t0,
                "find": t2-t1,
                "chains": t3-t2,
            },
            "counts": {
                "simple": simple,
                "super": super_,
                "insertion": insertion,
            },
            "bubbles": bubbles,
        }


    def memory(self,gfa_file):

        result=subprocess.run(
            [sys.executable,"-c",MEMORY_SCRIPT,gfa_file],
            capture_output=True,
            text=True
        )

        if result.returncode != 0:
            raise RuntimeError(f"BubbleGun memory run failed on {gfa_file}:\n{result.stderr}")

        baseline,peak=map(int,result.stdout.split())

        return {"peak": peak, "baseline": baseline}


    def close(self):

        logging.getLogger().removeHandler(self.warnings)
