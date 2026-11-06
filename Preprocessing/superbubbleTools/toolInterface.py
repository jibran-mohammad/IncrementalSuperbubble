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

    time      seconds of one run of the tool on gfa_file (the work a
              static tool must redo after every update; tool start-up
              should not be included)
    phases    {phase: seconds}, optional breakdown of the run
    counts    {bubble_type: number}, every type the tool reports
    bubbles   list of (bubble_type, end1, end2) with GFA segment
              names; filled only when keep_bubbles is True

    keep_bubbles must not change what "time" measures.

memory(gfa_file)   (optional; return None if not supported)

    Peak memory of the tool on gfa_file, measured in a separate process
    so the caller's own memory is not included. Return a dict:

    peak      peak resident memory (KB) of that process, comparable with
              PeakMemory(KB) of the C++ algorithms
    baseline  resident memory (KB) of the process before the graph is
              read (interpreter + tool code); peak - baseline is what
              the graph and its processing cost
"""


class SuperbubbleTool:

    name = None

    def run(self,gfa_file,keep_bubbles=False):

        raise NotImplementedError


    def memory(self,gfa_file):

        return None
