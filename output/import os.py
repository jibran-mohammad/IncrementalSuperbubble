import os
import re
import matplotlib.pyplot as plt
import pandas as pd
from glob import glob

# Config
input_files = glob("/home/dell/Research/IncSuperbubble/output/Human/MiniGraph/timeFile/GRCh38-464.gfa/*.txt")  # Folder with output files
k = 5  # Averaging window size

# Output containers
summary_data = []
plot_data = {}

# Helper function to parse file
def parse_file(filepath):
    with open(filepath) as f:
        lines = f.readlines()

    num_edges, num_nodes = map(int, lines[0].split())
    edge_data = []

    # Extract time and superbubbles
    for line in lines[1:-1]:
        parts = line.strip().split()
        if len(parts) < 3:
            continue
        edge_no = int(parts[0])
        time = float(parts[1])
        bubbles = int(parts[2])
        edge_data.append((edge_no, time, bubbles))

    # Peak memory is in the last line
    mem_line = lines[-1]
    match = re.search(r"PeakMemory\(KB\):\s*(\d+)", mem_line)
    peak_memory_kb = int(match.group(1)) if match else 0
    peak_memory_mb = peak_memory_kb / 1024  # Convert to MB

    total_time = sum([x[1] for x in edge_data])
    return total_time, peak_memory_mb, edge_data
    
    

# Process each file
for filepath in input_files:
    algo_name = os.path.basename(filepath).replace(".txt", "")
    total_time, peak_memory, edge_data = parse_file(filepath)

    summary_data.append([algo_name, total_time, peak_memory])
    plot_data[algo_name] = edge_data

# Create table
df_summary = pd.DataFrame(summary_data, columns=["Algorithm", "Total Time", "Peak Memory (MB)"])
latex_table = df_summary.to_latex(index=False, float_format="%.5f")
with open("summary_table.tex", "w") as f:
    f.write(latex_table)
print("\nLaTeX table written to 'summary_table.tex'")

# Plotting
plt.figure(figsize=(10, 6))

for algo, data in plot_data.items():
    # Smooth using window of k points
    smoothed_edges = []
    smoothed_times = []
    smoothed_bubbles = []

    for i in range(0, len(data), k):
        chunk = data[i:i+k]
        if not chunk:
            continue
        avg_edge = sum([x[0] for x in chunk]) / len(chunk)
        avg_time = sum([x[1] for x in chunk]) / len(chunk)
        avg_bubbles = sum([x[2] for x in chunk]) / len(chunk)
        smoothed_edges.append(avg_edge)
        smoothed_times.append(avg_time)
        smoothed_bubbles.append(avg_bubbles)

    plt.plot(smoothed_edges, smoothed_times, label=f"{algo}")
    for x, y, bubble in zip(smoothed_edges, smoothed_times, smoothed_bubbles):
        plt.text(x, y, str(int(bubble)), fontsize=8, ha='center', va='bottom', alpha=0.6)

plt.xlabel("Number of Edges")
plt.ylabel("Time")
plt.title("Time vs Edges with Superbubble Count")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("performance_plot.png")
plt.show()
