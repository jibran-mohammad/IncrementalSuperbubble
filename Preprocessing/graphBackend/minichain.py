import os
import shlex
from .base import GraphBackend
from ..utils import run_command

class MinichainBackend(GraphBackend):
    def __init__(self, threads, minichain_bin):
        self.threads = threads
        self.minichain_bin = minichain_bin

    def build_graph(self, ref_fasta, fasta_list, output_gfa, label=None):
        """
        Build a pangenome graph using minichain.

        label: optional string for logging (e.g. chromosome_1)
        """

        cmd = f"{shlex.quote(self.minichain_bin)} -t{self.threads} -cxggs {shlex.quote(ref_fasta)}"
        for fasta in fasta_list:
            if os.path.abspath(fasta) == os.path.abspath(ref_fasta):
                continue
            cmd += " " + shlex.quote(fasta)
        cmd += " > " + shlex.quote(output_gfa)
        
        print(f"Generating graph from {len(fasta_list)} FASTA files")
        run_command(cmd, label=label)
