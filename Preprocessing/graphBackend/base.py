from abc import ABC, abstractmethod

class GraphBackend(ABC):
    def __init__(self, threads):
        self.threads = threads

    @abstractmethod
    def build_graph(self, ref_fasta, fasta_list, output_gfa):
        """
        Must generate a GFA graph at output_gfa.
        """
        pass
