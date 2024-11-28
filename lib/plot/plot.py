import numpy as np
import matplotlib.pyplot as plt

# Function to load data and extract m and n
def load_data(filename):
    with open(filename, 'r') as file:
        # Read the first line for m and n
        first_line = file.readline().strip()
        n, m = map(int, first_line.split())
    # Load the rest of the file (edge, time, case), only taking edge and time columns
    edges, times = np.loadtxt(filename, skiprows=1, usecols=(0, 1), unpack=True)
    return m, n, edges, times

# Load data from each file with the new format
m1, n1, edges1, time1 = load_data('/home/dell/Research/IncSuperbubble-master/output/Email/SuperbubbleSimpleTimeOutput.txt')
m2, n2, edges2, time2 = load_data('/home/dell/Research/IncSuperbubble-master/output/Email/SuperbubbleincrementalTimeOutput.txt')
m3, n3, edges3, time3 = load_data('/home/dell/Research/IncSuperbubble-master/output/Email/SuperbubbleincrementalAVLTimeOutput.txt')

# Downsample real data more aggressively for smoother plots
downsample_rate = 50  # Use every 50th point for better smoothing
edges1_smooth, time1_smooth = edges1[::downsample_rate], time1[::downsample_rate]
edges2_smooth, time2_smooth = edges2[::downsample_rate], time2[::downsample_rate]
edges3_smooth, time3_smooth = edges3[::downsample_rate], time3[::downsample_rate]

# Plotting the actual data without complexity curves
plt.figure(figsize=(10, 6))
plt.plot(edges1_smooth, time1_smooth, label="File 1: O(m^2)", linewidth=0.5, alpha=0.7, color='red')
plt.plot(edges2_smooth, time2_smooth, label="File 2: O(m^2 * log(n))", linewidth=0.5, alpha=0.7, color='green')
plt.plot(edges3_smooth, time3_smooth, label="File 3: O(m * n * log(n))", linewidth=0.5, alpha=0.7, color='blue')

# Labeling the graph
plt.xlabel("Number of Edges (m)")
plt.ylabel("Time (seconds)")
plt.title("Time Complexity Comparisons with Smoothed Data Points")
plt.yscale('log')
plt.ylim([1e-5, 1])  # Set y-axis limit to focus on data range
plt.legend()
plt.grid(True)

# Show the plot
plt.show()
