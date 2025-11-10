import os
import sys

def parse_gfa(file_path):
    """
    Parses a GFA file and returns sets of oriented nodes and edges,
    and a dictionary of segment sequences.
    
    Nodes are represented as tuples: (segment_name, orientation).
    Edges are represented as tuples: (from_oriented_node, to_oriented_node).
    """
    nodes = set()
    edges = set()
    sequences = {}

    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            line_type = parts[0]

            if line_type == 'S':
                # Add both forward and reverse oriented nodes
                segment_name = parts[1]
                sequence = parts[2]
                
                nodes.add((segment_name, '+'))
                nodes.add((segment_name, '-'))
                sequences[segment_name] = sequence

            elif line_type == 'L':
                if len(parts) >= 5:
                    from_seg, from_or = parts[1], parts[2]
                    to_seg, to_or = parts[3], parts[4]

                    # Add the oriented nodes involved in the link
                    nodes.add((from_seg, from_or))
                    nodes.add((to_seg, to_or))

                    # Add the directed edge between oriented nodes
                    edges.add(((from_seg, from_or), (to_seg, to_or)))
                    
    return nodes, edges, sequences

def analyze_graph_differences(graph1_path, graph2_path):
    """
    Analyzes and reports the differences between two GFA graphs.
    """
    print(f"Parsing '{graph1_path}'...")
    nodes1, edges1, seqs1 = parse_gfa(graph1_path)
    
    print(f"Parsing '{graph2_path}'...")
    nodes2, edges2, seqs2 = parse_gfa(graph2_path)

    print("\n--- Graph Difference Analysis ---")
    print(f"Graph 1: {os.path.basename(graph1_path)}")
    print(f"Graph 2: {os.path.basename(graph2_path)}")
    print("-" * 35)

    # 1. Sequence Differences
    print("\n1. Segments with Sequence Mismatches:")
    mismatched_segments = []
    common_segments = set(seqs1.keys()) & set(seqs2.keys())
    for segment_name in common_segments:
        if seqs1[segment_name] != seqs2[segment_name]:
            mismatched_segments.append(segment_name)
    
    if mismatched_segments:
        print(f"   Found {len(mismatched_segments)} segments with different sequences:")
        for name in sorted(mismatched_segments):
            print(f"   - {name}")
    else:
        print("   No sequence mismatches found for common segments.")

    # 2. Node (Oriented Segment) Differences
    print("\n2. Unique Nodes (Oriented Segments):")
    nodes_only_in_g1 = nodes1 - nodes2
    nodes_only_in_g2 = nodes2 - nodes1
    
    print(f"   Nodes in Graph 1 only ({len(nodes_only_in_g1)}):")
    if nodes_only_in_g1:
        for node in sorted(list(nodes_only_in_g1)):
            print(f"   - Segment: {node[0]}, Orientation: {node[1]}")
    
    print(f"\n   Nodes in Graph 2 only ({len(nodes_only_in_g2)}):")
    if nodes_only_in_g2:
        for node in sorted(list(nodes_only_in_g2)):
            print(f"   - Segment: {node[0]}, Orientation: {node[1]}")
    
    # 3. Edge Differences
    print("\n3. Unique Edges:")
    edges_only_in_g1 = edges1 - edges2
    edges_only_in_g2 = edges2 - edges1
    
    print(f"   Edges in Graph 1 only ({len(edges_only_in_g1)}):")
    if edges_only_in_g1:
        for edge in sorted(list(edges_only_in_g1)):
            from_node, to_node = edge
            print(f"   - From ({from_node[0]}{from_node[1]}) -> To ({to_node[0]}{to_node[1]})")
    
    print(f"\n   Edges in Graph 2 only ({len(edges_only_in_g2)}):")
    if edges_only_in_g2:
        for edge in sorted(list(edges_only_in_g2)):
            from_node, to_node = edge
            print(f"   - From ({from_node[0]}{from_node[1]}) -> To ({to_node[0]}{to_node[1]})")

if __name__ == "__main__":
    # Example usage: Replace with your actual file paths
    graph1_path = "/home/dell/Research/Superbubble/Experiment/IncSuperbubble/output/graphs_incremental/minigraph/Incremental/graph1.gfa"
    graph2_path = "/home/dell/Research/Superbubble/Experiment/IncSuperbubble/output/graphs_incremental/minigraph/Incremental/graph2.gfa"

    if not os.path.exists(graph1_path):
        sys.exit(f"Error: Graph 1 file not found at '{graph1_path}'")
    if not os.path.exists(graph2_path):
        sys.exit(f"Error: Graph 2 file not found at '{graph2_path}'")

    analyze_graph_differences(graph1_path, graph2_path)