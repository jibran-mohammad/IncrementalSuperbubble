import pandas as pd
import matplotlib
# Set the backend to 'Agg' BEFORE importing pyplot
# This is a non-interactive backend, perfect for saving files
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import sys
import os # Import the os module to handle paths

def plot_haplotype_data(filename, output_dir):
    """
    Reads the haplotype comparison CSV and generates three plots,
    saving them to the specified output_dir.
    
    Plots generated:
    1. Runtimes (log scale) vs. Haplotype Step
    2. Runtimes (log scale) vs. Number of Edges
    3. Number of Superbubbles vs. Number of Edges
    """
    
    # --- File Reading ---
    try:
        df = pd.read_csv(filename)
        print(f"Successfully loaded '{filename}'. Found {len(df)} rows.")
    except FileNotFoundError:
        print(f"Error: The file '{filename}' was not found.")
        return
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return

    # --- Create output directory if it doesn't exist ---
    try:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            print(f"Created output directory: {output_dir}")
    except Exception as e:
        print(f"Error: Could not create output directory '{output_dir}': {e}")
        return

    # --- Data Cleaning ---
    # The 'Step' column contains 'Base', which we'll replace with '0' for plotting.
    if 'Step' in df.columns:
        df['Step'] = df['Step'].replace('Base', 0)
    else:
        print("Warning: 'Step' column not found.")
        return

    # Define columns that should be numeric
    cols_to_numeric = [
        'Step', 'Nodes', 'Edges', 'Superbubbles', 
        'BST_INC_DFS_Time_s', 'INC_DFS_Time_s', 'Simple_Time_s'
    ]
    
    # Ensure all necessary columns exist before trying to convert
    missing_cols = [col for col in cols_to_numeric if col not in df.columns]
    if missing_cols:
        print(f"Error: The CSV is missing the following required columns: {missing_cols}")
        return

    # Convert all plotting columns to numeric, coercing errors to NaN
    for col in cols_to_numeric:
        df[col] = pd.to_numeric(df[col], errors='coerce')
        
    # Drop any rows that had conversion errors
    original_rows = len(df)
    df = df.dropna(subset=cols_to_numeric)
    cleaned_rows = len(df)
    if original_rows > cleaned_rows:
        print(f"Data cleaning: Removed {original_rows - cleaned_rows} rows with invalid data.")
    print("Data cleaning complete.")

    # --- Plot 1: Runtimes vs. Haplotype Step ---
    print("Generating Plot 1: Runtimes vs. Haplotype Step")
    df_step_sorted = df.sort_values(by='Step')
    
    plt.figure(figsize=(12, 7))
    plt.plot(df_step_sorted['Step'], df_step_sorted['BST_INC_DFS_Time_s'], marker='o', linestyle='-', label='AVLDFS')
    plt.plot(df_step_sorted['Step'], df_step_sorted['INC_DFS_Time_s'], marker='s', linestyle='--', label='DFS')
    plt.plot(df_step_sorted['Step'], df_step_sorted['Simple_Time_s'], marker='^', linestyle=':', label='Simple')
    
    plt.title('Algorithm Runtimes vs. Haplotype Step', fontsize=16)
    plt.xlabel('Haplotype Step', fontsize=12)
    plt.ylabel('Time (s) [Log Scale]', fontsize=12)
    plt.yscale('log') # Set y-axis to log scale
    plt.legend()
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tight_layout()
    
    # Save the figure instead of showing it
    plot1_filename = os.path.join(output_dir, "plot1_runtimes_vs_step.png")
    plt.savefig(plot1_filename)
    print(f"Saved Plot 1 to {plot1_filename}")
    plt.close() # Close the figure to free memory

    # --- Plot 2: Runtimes vs. Number of Edges ---
    print("Generating Plot 2: Runtimes vs. Number of Edges")
    df_edges_sorted = df.sort_values(by='Edges')
    
    plt.figure(figsize=(12, 7))
    plt.plot(df_edges_sorted['Edges'], df_edges_sorted['BST_INC_DFS_Time_s'], marker='o', linestyle='-', label='AVLDFS')
    plt.plot(df_edges_sorted['Edges'], df_edges_sorted['INC_DFS_Time_s'], marker='s', linestyle='--', label='DFS')
    plt.plot(df_edges_sorted['Edges'], df_edges_sorted['Simple_Time_s'], marker='^', linestyle=':', label='Simple')
    
    plt.title('Algorithm Runtimes vs. Number of Edges', fontsize=16)
    plt.xlabel('Number of Edges', fontsize=12)
    plt.ylabel('Time (s) [Log Scale]', fontsize=12)
    plt.yscale('log') # Set y-axis to log scale
    plt.legend()
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tight_layout()
    
    # Save the figure instead of showing it
    plot2_filename = os.path.join(output_dir, "plot2_runtimes_vs_edges.png")
    plt.savefig(plot2_filename)
    print(f"Saved Plot 2 to {plot2_filename}")
    plt.close() # Close the figure

    # --- Plot 3: Number of Superbubbles vs. Number of Edges ---
    print("Generating Plot 3: Number of Superbubbles vs. Number of Edges")
    # We can reuse df_edges_sorted from Plot 2

    plt.figure(figsize=(12, 7))
    # This plot shows the relationship between two graph metrics
    plt.plot(df_edges_sorted['Edges'], df_edges_sorted['Superbubbles'], marker='o', linestyle='-', label='Superbubbles Count', color='b')

    plt.title('Graph Properties: Number of Superbubbles vs. Number of Edges', fontsize=16)
    plt.xlabel('Number of Edges', fontsize=12)
    plt.ylabel('Number of Superbubbles', fontsize=12)
    # Using a linear scale here to show the direct relationship
    plt.legend()
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tight_layout()
    
    # Save the figure instead of showing it
    plot3_filename = os.path.join(output_dir, "plot3_superbubbles_vs_edges.png")
    plt.savefig(plot3_filename)
    print(f"Saved Plot 3 to {plot3_filename}")
    plt.close() # Close the figure

    print(f"All plots generated and saved to {output_dir}.")

if __name__ == "__main__":
    """
    Main entry point for the script.
    Expects two command-line arguments:
    1. The path to the input CSV file.
    2. The path to the output directory for plots.
    
    Example:
    python plot_haplotype_data.py data/my_data.csv plots/
    """
    
    # Ensure the user has the required packages
    # You can install them with:
    # pip install pandas matplotlib
    
    # Check if a filename AND output dir were provided as command-line arguments
    if len(sys.argv) > 2:
        # sys.argv[0] is the script name, sys.argv[1] is the input file, sys.argv[2] is the output dir
        filepath = sys.argv[1]
        output_path = sys.argv[2]
        plot_haplotype_data(filepath, output_path)
    else:
        # If not enough arguments are provided, show the user how to run the script
        print("Usage: python plot_haplotype_data.py <path_to_your_csv_file> <output_directory>")
        print("Example: python plot_haplotype_data.py data/my_data.csv plots/")
        print("Please provide both the input file path and the desired output directory.")
        sys.exit(1) # Exit with an error code