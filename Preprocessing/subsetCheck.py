import os
import sys

def parse_gfa_by_sequence(file_path):
    """
    Parses a GFA file and returns sets of unique sequences and sequence-based edges.
    
    This function ignores segment names and treats the DNA sequence itself as the
    identifier for a node, to check for a sequence-based subset relationship.
    
    Returns:
        tuple: (sequences_set, edges_set)
    """
    sequences = {}
    oriented_sequence_to_links = {}

    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            line_type = parts[0]

            if line_type == 'S':
                segment_name = parts[1]
                sequence = parts[2]
                sequences[segment_name] = sequence

            elif line_type == 'L':
                if len(parts) >= 5:
                    from_seg, from_or = parts[1], parts[2]
                    to_seg, to_or = parts[3], parts[4]

                    # Store links based on oriented segment names first
                    if (from_seg, from_or) not in oriented_sequence_to_links:
                        oriented_sequence_to_links[(from_seg, from_or)] = set()
                    oriented_sequence_to_links[(from_seg, from_or)].add((to_seg, to_or))

    # Now convert the links from segment names to sequences
    unique_sequences = set(sequences.values())
    sequence_based_edges = set()

    for (from_seg, from_or), to_nodes in oriented_sequence_to_links.items():
        if from_seg in sequences:
            from_seq = sequences[from_seg]
            for to_seg, to_or in to_nodes:
                if to_seg in sequences:
                    to_seq = sequences[to_seg]
                    sequence_based_edges.add(((from_seq, from_or), (to_seq, to_or)))

    return unique_sequences, sequence_based_edges

def compare_graphs_by_sequence(graph1_path, graph2_path):
    """
    Compares two graphs based on their unique sequences and sequence-based edges.
    """
    print(f"Parsing {graph1_path} by sequence...")
    seqs1, edges1 = parse_gfa_by_sequence(graph1_path)
    
    print(f"Parsing {graph2_path} by sequence...")
    seqs2, edges2 = parse_gfa_by_sequence(graph2_path)
    
    # Check for sequence-based subset relationship
    is_sequence_subset = seqs1.issubset(seqs2)
    is_edge_subset = edges1.issubset(edges2)

    if is_sequence_subset and is_edge_subset:
        print(f"\n'{graph1_path}' IS a sequence-based subset of '{graph2_path}'.")
        
        # Find new sequences and edges
        new_sequences = seqs2 - seqs1
        new_edges = edges2 - edges1
        
        print("\n--- Summary of Sequence-based Changes ---")
        print(f"Total unique sequences in Graph 1: {len(seqs1)}")
        print(f"Total unique sequences in Graph 2: {len(seqs2)}")
        print(f"New sequences added: {len(new_sequences)}")
        print(f"\nTotal edges in Graph 1: {len(edges1)}")
        print(f"Total edges in Graph 2: {len(edges2)}")
        print(f"New edges added: {len(new_edges)}")
    else:
        print(f"\n'{graph1_path}' is NOT a sequence-based subset of '{graph2_path}'.")
        if not is_sequence_subset:
            print(f"{len(seqs1 - seqs2)} unique sequences from Graph 1 not found in Graph 2.")
        if not is_edge_subset:
            print(f"{len(edges1 - edges2)} edges from Graph 1 not found in Graph 2.")

if __name__ == "__main__":
    # Example usage with your file paths
    # Replace these with the actual paths to your graph files
    graph1_path = "/home/dell/Research/Superbubble/Experiment/IncSuperbubble/output/graphs_incremental/minichain/Sanitized/graph1_sanitized.gfa"
    graph2_path = "/home/dell/Research/Superbubble/Experiment/IncSuperbubble/output/graphs_incremental/minichain/Sanitized/graph2_sanitized.gfa"

    if not os.path.exists(graph1_path):
        sys.exit(f"Error: Graph 1 file not found at {graph1_path}")
    if not os.path.exists(graph2_path):
        sys.exit(f"Error: Graph 2 file not found at {graph2_path}")

    compare_graphs_by_sequence(graph1_path, graph2_path)