import matplotlib.pyplot as plt

# File path (update to your file path)
file_path = "../../output/Yeast/timeFile/SuperbubbleincrementalDFSAVLTreeOutput.txt"

# Read data from the file
with open(file_path, 'r') as file:
    lines = file.readlines()

# Extract the number of edges and nodes from the first line
num_edges, num_nodes = map(int, lines[0].strip().split())

# Extract data for plotting
edge_numbers = []
superbubble_counts = []

for line in lines[1:]:
    parts = line.strip().split()
    if len(parts) >= 3:
        edge_numbers.append(int(parts[0]))  # Edge number
        superbubble_counts.append(int(parts[2]))  # Superbubble count

# Plot the graph
plt.plot(edge_numbers, superbubble_counts, marker='o', linestyle='-', color='b', linewidth=0.1, label='Superbubbles')

# Add labels, title, and legend
plt.xlabel("Number of Edges")
plt.ylabel("Number of Superbubbles")
plt.title(f"Superbubble Growth with Edge Insertions (Nodes: {num_nodes})")
plt.legend()

# Annotate the number of nodes on the plot
plt.text(0.05, 0.95, f"Number of Nodes: {num_nodes}", transform=plt.gca().transAxes,
         fontsize=10, bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray'))

# Show the grid and plot
plt.grid(True)
plt.tight_layout()
plt.show()
