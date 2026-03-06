"""
toolInterface.py

Abstract interface for superbubble tools.

Any tool must implement:

run(gfa_file)

Return:

runtime
superbubble_count
"""


class SuperbubbleTool:

    def run(self,gfa_file):

        raise NotImplementedError