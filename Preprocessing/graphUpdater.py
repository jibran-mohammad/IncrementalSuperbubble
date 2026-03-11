"""
graphUpdater.py

Maintains an evolving GFA graph in memory.

Used for:

- perEdge updates
- perBatch updates

The base graph is read once.
Edges are added incrementally.

Temporary GFA files are created
for running superbubble tools.
"""

import tempfile


class GraphUpdater:

    def __init__(self,base_gfa):

        with open(base_gfa) as f:
            self.base_lines=f.readlines()

        self.current_edges=set()


    def add_edge(self,edge):

        if edge in self.current_edges:
            return

        self.current_edges.add(edge)


    def write_temp_gfa(self):

        temp_name = "/tmp/temp_graph.gfa"

        with open(temp_name,'w') as f:

            for line in self.base_lines:
                f.write(line)

            for edge in self.current_edges:

                u,v=edge.split("->")

                u=u[:-1]
                v=v[:-1]

                f.write(f"L\t{u}\t+\t{v}\t+\t0M\n")

        return temp_name