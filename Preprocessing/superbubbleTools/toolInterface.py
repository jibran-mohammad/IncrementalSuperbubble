"""
toolInterface.py

Abstract interface for superbubble tools.

Any tool must define:

name
    Short lowercase tool name (e.g. "bubblegun"). Used as the
    output directory and file-name prefix, so results land in
    output/<dataset>/<name>/<mode>/.

run(gfa_file, keep_bubbles=False)

    Return a dict:

    time      wall-clock seconds of one run of the tool on gfa_file
    counts    {bubble_type: number}, every type the tool reports
    bubbles   list of (bubble_type, end1, end2) with GFA segment
              names; filled only when keep_bubbles is True

    keep_bubbles must not change what "time" measures.
"""


class SuperbubbleTool:

    name = None

    def run(self,gfa_file,keep_bubbles=False):

        raise NotImplementedError
