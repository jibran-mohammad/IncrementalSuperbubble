import numpy as np
import matplotlib.pyplot as plt

# Function to load data and extract m and n
def load_data(filename):
    with open(filename, 'r') as file:
        # Read the first line for m and n
        first_line = file.readline().strip()
        n, m = map(int, first_line.split())
    # Load the rest of the file (edge, time), taking the columns for edges and times
    edges, times = np.loadtxt(filename, skiprows=1, usecols=(0, 1), unpack=True)
    return m, n, edges, times

# Function to compute average in chunks
def chunk_average(x, y, chunk_size):
    averaged_x = []
    averaged_y = []
    for i in range(0, len(x), chunk_size):
        chunk_x = x[i:i + chunk_size]
        chunk_y = y[i:i + chunk_size]
        if len(chunk_x) > 0:  # Avoid empty chunks
            averaged_x.append(np.mean(chunk_x))
            averaged_y.append(np.mean(chunk_y))
    return np.array(averaged_x), np.array(averaged_y)

# Chunk size for averaging
chunk_size = 150  # Adjust as needed

# Load data from files
m1, n1, edges1, time1 = load_data('/home/dell/Research/IncSuperbubble/output/Yeast/timeFile/SuperbubbleSimpleTimeOutput.txt')
m2, n2, edges2, time2 = load_data('/home/dell/Research/IncSuperbubble/output/Yeast/timeFile/SuperbubbleincrementalDFSOutput.txt')
m3, n3, edges3, time3 = load_data('/home/dell/Research/IncSuperbubble/output/Yeast/timeFile/SuperbubbleincrementalDFSAVLTreeOutput.txt')
m4, n4, edges4, time4 = load_data('/home/dell/Research/IncSuperbubble/output/Yeast/timeFile/SuperbubbleincrementalH1output.txt')

# Calculate averages in chunks
edges1_avg, time1_avg = chunk_average(edges1, time1, chunk_size)
edges2_avg, time2_avg = chunk_average(edges2, time2, chunk_size)
edges3_avg, time3_avg = chunk_average(edges3, time3, chunk_size)
edges4_avg, time4_avg = chunk_average(edges4, time4, chunk_size)

# Plotting the data with averaged points
plt.figure(figsize=(10, 6))
plt.plot(edges1_avg, time1_avg, label=f"File 1: O(m^2), n={m1}", linewidth=1.0, alpha=0.7, color='red')
plt.plot(edges2_avg, time2_avg, label=f"File 2: O(m^2 * log(n)), n={m2}", linewidth=1.0, alpha=0.7, color='green')
plt.plot(edges3_avg, time3_avg, label=f"File 3: O(m * n * log(n)), n={m3}", linewidth=1.0, alpha=0.7, color='blue')
plt.plot(edges4_avg, time4_avg, label=f"File 4: O(m * n * log(n) H1), n={m4}", linewidth=1.0, alpha=0.7, color='black')


# Labeling the graph
plt.xlabel("Number of Edges (m)")
plt.ylabel("Time (seconds)")
plt.title("Time Complexity Comparisons with Averaged Data Points")
plt.yscale('log')
plt.ylim([1e-5, 1])  # Set y-axis limit to focus on data range
plt.legend()
plt.grid(True)

# Show the plot
plt.show()
