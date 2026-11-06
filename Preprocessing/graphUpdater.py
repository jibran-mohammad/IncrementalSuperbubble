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

Nodes that are not in the base graph get an S line with
sequence '*' when an edge first uses them; GFA readers such
as BubbleGun skip L lines whose segments have no S line.
"""

import os
import tempfile


class GraphUpdater:

    def __init__(self,base_gfa):

        with open(base_gfa) as f:
            self.base_lines=f.readlines()

        self.base_nodes={
            line.split("\t")[1]
            for line in self.base_lines
            if line.startswith("S\t")
        }

        # dicts used as insertion-ordered sets: the temp GFA lists nodes
        # and edges in the order they were added, so it is identical on
        # every run (set order changes between Python processes, and
        # BubbleGun's choice of bubble source/sink depends on line order)
        self.current_edges={}
        self.new_nodes={}

        # Private temp file, so parallel runs do not overwrite each other.
        fd,self.temp_name=tempfile.mkstemp(prefix="temp_graph_",suffix=".gfa")
        os.close(fd)


    def add_edge(self,edge):

        if edge in self.current_edges:
            return

        self.current_edges[edge]=None

        u,v=edge.split("->")

        for node in (u[:-1],v[:-1]):

            if node not in self.base_nodes:
                self.new_nodes[node]=None


    def write_temp_gfa(self):

        temp_name = self.temp_name

        with open(temp_name,'w') as f:

            for line in self.base_lines:
                f.write(line)

            for node in self.new_nodes:
                f.write(f"S\t{node}\t*\n")

            for edge in self.current_edges:

                u,v=edge.split("->")

                u=u[:-1]
                v=v[:-1]

                f.write(f"L\t{u}\t+\t{v}\t+\t0M\n")

        return temp_name


    def cleanup(self):

        if os.path.exists(self.temp_name):
            os.remove(self.temp_name)
