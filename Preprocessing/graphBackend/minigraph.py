import os
import shlex
from .base import GraphBackend
from ..utils import run_command

class MinigraphBackend(GraphBackend):
    def __init__(self, threads, minigraph_bin):
        self.threads = threads
        self.minigraph_bin = minigraph_bin

    def build_graph(self, ref_fasta, fasta_list, output_gfa, label=None):
        # minigraph usually wants reference first
        cmd = f"{shlex.quote(self.minigraph_bin)} -t {self.threads} --inv=no -cxggs {shlex.quote(ref_fasta)}"
        for fasta in fasta_list:
            if os.path.abspath(fasta) == os.path.abspath(ref_fasta):
                continue
            cmd += " " + shlex.quote(fasta)
        cmd += f" > {shlex.quote(output_gfa)}"

        print(f"Generating graph from {len(fasta_list)} FASTA files")
        run_command(cmd, label=label)
