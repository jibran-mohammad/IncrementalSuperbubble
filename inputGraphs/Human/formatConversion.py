def convert_gfa_to_dag_format(gfa_file, output_file="converted_dag.dag"):
    """
    Convert GFA graph file to the DAG format expected by your algorithm.
    
    """
    edges = []  # List to store edges (source, target)
    all_nodes = set()  # Set to track all unique nodes

    # Read the GFA file and extract edges
    with open(gfa_file, 'r') as file:
        for line in file:
            if line.startswith("L"):  # Only process Link lines
                parts = line.strip().split('\t')
                source = int(parts[1])  # Extract Source node
                target = int(parts[3])  # Extract target node
                edges.append((source, target))
                all_nodes.update([source, target])  # Add nodes to the set

    # Find the highest node number and determine the new node number
    highest_node = max(all_nodes)
    print(highest_node)
    replacement_node = highest_node + 1
    print(replacement_node)

    # Replace node 1 with the new node number in the edges
    updated_edges = [(replacement_node if src == 1 else src, replacement_node if tgt == 1 else tgt) 
        for src, tgt in edges
    ]            
    
    # Write the edges to a new output file in the required format
    with open(output_file, 'w') as new_file:
        new_file.write(f"{len(edges)}\n")  # First line: Number of edges
        for edge in updated_edges:
            new_file.write(f"{edge[0]} {edge[1]}\n")  # Write each edge
    
    print(f"DAG format successfully written to {output_file}")


convert_gfa_to_dag_format("minigraphcactus/hprc-v1.1-mc-grch38.gfa", "minigraphcactus/hprc-v1.1-mc-grch38Converted.gfa")
