import os
import sys
import glob

def find_duplicate_sequences(gfa_file_path):
    """
    Scans a GFA file and returns a dictionary of sequences that appear
    with multiple different segment names.
    
    Returns:
        dict: A dictionary mapping a duplicate sequence to a list of
              segment names that share it.
    """
    sequences = {}
    duplicates = {}

    try:
        with open(gfa_file_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) > 2 and parts[0] == 'S':
                    segment_name = parts[1]
                    sequence = parts[2]
                    
                    if sequence in sequences and sequences[sequence] != segment_name:
                        if sequence not in duplicates:
                            duplicates[sequence] = [sequences[sequence]]
                        if segment_name not in duplicates[sequence]:
                            duplicates[sequence].append(segment_name)
                    else:
                        sequences[sequence] = segment_name
                        
    except FileNotFoundError:
        print(f"Error: File not found at '{gfa_file_path}'")
        return None
        
    return duplicates

def count_nodes_and_edges_with_duplicate_sequences(gfa_file_path, duplicates):
    """
    Counts the total number of oriented nodes and edges associated with
    segments that have a duplicate sequence.
    
    Args:
        gfa_file_path (str): Path to the GFA file.
        duplicates (dict): Dictionary from find_duplicate_sequences.
    
    Returns:
        tuple: (total_nodes, total_edges)
    """
    duplicate_segment_names = set()
    for names in duplicates.values():
        duplicate_segment_names.update(names)

    node_count = 0
    edge_count = 0
    
    try:
        with open(gfa_file_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                line_type = parts[0]

                if line_type == 'S':
                    segment_name = parts[1]
                    if segment_name in duplicate_segment_names:
                        node_count += 2  # One for '+' and one for '-' orientation
                
                elif line_type == 'L':
                    if len(parts) >= 5:
                        from_seg, to_seg = parts[1], parts[3]
                        if from_seg in duplicate_segment_names or to_seg in duplicate_segment_names:
                            edge_count += 1
                            
    except FileNotFoundError:
        print(f"Error: File not found at '{gfa_file_path}'")
        return None
    
    return node_count, edge_count


if __name__ == "__main__":
    # Example usage: Replace with your actual file path
    gfa_file = "/home/dell/Research/Superbubble/Experiment/IncSuperbubble/output/graphs_incremental/minichain/Incremental/graph1.gfa"

    if not os.path.exists(gfa_file):
        sys.exit(f"Error: GFA file not found at '{gfa_file}'")

    print(f"Scanning '{os.path.basename(gfa_file)}' for duplicate sequences...")
    duplicates = find_duplicate_sequences(gfa_file)

    if not duplicates:
        print("No duplicate sequences found.")
    else:
        print("\n--- Found duplicate sequences with different segment names ---")
        for seq, names in duplicates.items():
            print(f"Sequence: '{seq}'")
            print(f"  Shared by segment names: {', '.join(names)}")
            print("-" * 20)
        
        nodes, edges = count_nodes_and_edges_with_duplicate_sequences(gfa_file, duplicates)
        
        print("\n--- Summary of duplicate-sequence-related elements ---")
        print(f"Total sequences with duplicates: {len(duplicates)}")
        print(f"Total oriented nodes with these sequences: {nodes}")
        print(f"Total edges associated with these nodes: {edges}")