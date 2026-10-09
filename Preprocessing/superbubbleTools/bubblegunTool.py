"""
bubblegunTool.py

BubbleGun implementation of the superbubble tool interface.

BubbleGun reports three bubble types (simple, super, insertion);
all of them are superbubbles, so all three are counted.

BubbleGun is run inside a private temporary directory so that its
log.log does not end up in the caller's working directory.
"""

import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

from .toolInterface import SuperbubbleTool


COUNT_PATTERNS = {
    "simple": r"number of Simple Bubbles is (\d+)",
    "super": r"number of Superbubbles is (\d+)",
    "insertion": r"number of insertions is (\d+)",
}


def find_bubblegun():
    """BubbleGun on PATH, else the one installed next to this Python."""

    exe = shutil.which("BubbleGun")

    if exe:
        return exe

    exe = os.path.join(os.path.dirname(sys.executable), "BubbleGun")

    if os.path.exists(exe):
        return exe

    raise FileNotFoundError(
        "BubbleGun executable not found (install it in the venv: pip install BubbleGun)"
    )


class BubbleGunTool(SuperbubbleTool):

    name = "bubblegun"

    def __init__(self):

        self.exe = find_bubblegun()
        self.work_dir = tempfile.mkdtemp(prefix="bubblegun_")


    def _bubblegun(self,gfa_file,extra_args=()):

        result=subprocess.run(
            [self.exe,"-g",os.path.abspath(gfa_file),"bchains",*extra_args],
            capture_output=True,
            text=True,
            cwd=self.work_dir
        )

        if result.returncode != 0:
            raise RuntimeError(
                f"BubbleGun failed on {gfa_file}:\n{result.stderr}"
            )

        return result


    def run(self,gfa_file,keep_bubbles=False):

        # Timed run: no JSON output, so the time is only BubbleGun's own work.
        start=time.time()

        result=self._bubblegun(gfa_file)

        end=time.time()

        counts={}

        for bubble_type,pattern in COUNT_PATTERNS.items():

            match=re.search(pattern,result.stdout)

            counts[bubble_type]=int(match.group(1)) if match else 0

        bubbles=[]

        # Separate, untimed run to get the bubble list.
        if keep_bubbles and sum(counts.values()) > 0:

            json_path=os.path.join(self.work_dir,"bubbles.json")

            self._bubblegun(gfa_file,("--bubble_json",json_path))

            with open(json_path) as f:
                chains=json.load(f)

            for chain in chains.values():
                for bubble in chain["bubbles"]:
                    end1,end2=bubble["ends"]
                    bubbles.append((bubble["type"],end1,end2))

            os.remove(json_path)

        return {
            "time": end-start,
            "counts": counts,
            "bubbles": bubbles,
        }


    def close(self):

        shutil.rmtree(self.work_dir,ignore_errors=True)
