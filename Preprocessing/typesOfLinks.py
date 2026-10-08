import os

def categorize_gfa_links(gfa_file):
    """
    Categorizes the links in a GFA file based on the orientation of the segments they connect.
    Args:
        gfa_file (str): The path to the input GFA file.
    Returns:
        A dictionary with counts for each link type.
    """
    link_counts = {
        '+ +': 0,
        '+ -': 0,
        '- +': 0,
        '- -': 0
    }
    
    with open(gfa_file, 'r') as f:
        for line in f:
            if line.startswith('L'):
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    from_or = parts[2]
                    to_or = parts[4]
                    
                    link_type = f"{from_or} {to_or}"
                    if link_type in link_counts:
                        link_counts[link_type] += 1
                        
    return link_counts

if __name__ == "__main__":
    # Example usage: Replace with your file path
    file_path = "/home/dell/Research/Superbubble/Experiment/IncSuperbubble/output/graphs_incremental/minigraph/Incremental/graph1.gfa"

    if not os.path.exists(file_path):
        print(f"Error: File not found at {file_path}")
    else:
        counts = categorize_gfa_links(file_path)
        
        total_links = sum(counts.values())
        if total_links == 0:
            print(f"No links found in {file_path}")
        else:
            print(f"Link categorization for {file_path}:")
            print("-" * 30)
            print(f"Total links: {total_links}")
            print("-" * 30)
            for link_type, count in counts.items():
                percentage = (count / total_links) * 100
                print(f"'{link_type}': {count} links ({percentage:.2f}%)")
            print("-" * 30)