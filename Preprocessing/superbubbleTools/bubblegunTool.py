"""
bubblegunTool.py

BubbleGun implementation of the superbubble tool interface.
"""

import subprocess
import time
import re

from .toolInterface import SuperbubbleTool


class BubbleGunTool(SuperbubbleTool):

    def run(self,gfa_file):

        start=time.time()

        result=subprocess.run(
            ["BubbleGun","-g",gfa_file,"bchains"],
            capture_output=True,
            text=True
        )

        end=time.time()

        runtime=end-start

        superbubbles=0

        match=re.search(
            r"Superbubbles is (\d+)",
            result.stdout
        )

        if match:
            superbubbles=int(match.group(1))

        return runtime,superbubbles