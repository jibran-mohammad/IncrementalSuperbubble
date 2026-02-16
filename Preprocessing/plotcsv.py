import pandas as pd
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import sys
import os
import numpy as np # Import numpy for bar chart positioning
from matplotlib.lines import Line2D
import argparse

# Define the components and their corresponding column names
# --- MODIFIED: 'Insert' component has been removed from this dictionary ---
TIME_COMPONENTS = {
    'Insert': ['SimpleInsertTime', 'INCInsertTime', 'BSTInsertTime'],
    'DFS': ['SimpleDFSTime', 'INCDFSTime', 'BSTDFSTime'],
    'MinI': ['SimpleMinITime', 'INCMinITime', 'BSTMinITime'],
    'MaxO': ['SimpleMaxOTime', 'INCMaxOTime', 'BSTMaxOTime'],
    'SBCompute': ['SimpleSBComputeTime', 'INCSBComputeTime', 'BSTSBComputeTime'],
    'Total': ['SimpleTotalTime', 'INCTotalTime', 'BSTTotalTime']
}

# Define color palette for the stacked bar chart components
# --- MODIFIED: 'Insert' color has been removed entirely ---
STACK_COLORS = {
    'Insert': '#4C72B0',
    'DFS': '#55A868',         # Green
    'MinI': '#C44E52',        # Red
    'MaxO': '#8172B2',        # Purple
    'SBCompute': '#CCB974',   # Yellow/Brown
}

# Define line styles for the three versions in line plots
VERSION_STYLES = {
    'Simple': {'color': 'red', 'marker': 'o', 'linestyle': '-',  "linewidth": 0.7, "markersize": 2},
    'INC': {'color': 'blue', 'marker': 's', 'linestyle': '--',  "linewidth": 0.7, "markersize": 2},
    'BST': {'color': 'green', 'marker': '^', 'linestyle': '-.',  "linewidth": 0.7, "markersize": 2},
}

def plot_total_time_bar(df, output_dir, dataset_name="dataset"):
    """
    Plot total time (sum across all steps excluding Base)
    for Simple, INC, BST as a single bar chart.
    """

    print("Generating Plot: Total Time Comparison Bar Chart...")

    df_plot = df[df["Step"] != 0].copy()

    versions = ["Simple", "INC", "BST"]
    totals = []

    for v in versions:
        col = f"{v}TotalTime"
        if col not in df_plot.columns:
            raise KeyError(f"Missing column: {col}")
        totals.append(df_plot[col].sum())

    fig, ax = plt.subplots(figsize=(6, 5))

    bars = ax.bar(versions, totals, color=["red", "blue", "green"],
                  edgecolor="black", linewidth=0.6)

    ax.set_ylabel("Total Time (ms)")
    ax.set_title(f"Total Time Comparison — {dataset_name}")
    ax.grid(axis="y", linestyle="--", alpha=0.6)

    # Add value labels above bars
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2,
                height * 1.01,
                f"{height:.2f}",
                ha='center', va='bottom', fontsize=9)

    plt.tight_layout()

    filename = os.path.join(output_dir, f"{dataset_name}_plot_total_time_bar.png")
    plt.savefig(filename, dpi=200, bbox_inches="tight")
    print(f"Saved total time bar plot to {filename}")
    plt.close(fig)


def clean_and_normalize_df(df):
    """
    1) Ensure Step is numeric and 'Base'->0
    2) Normalize expected column names (Simple/INC/BST + component)
    3) Coerce time columns to numeric and fill NaN with 0
    4) Fail fast with a readable error if a required column is missing
    """
    # 1) normalize Step
    if 'Step' in df.columns:
        df['Step'] = df['Step'].replace({'Base': 0})
        df['Step'] = pd.to_numeric(df['Step'], errors='coerce').fillna(0).astype(int)
    else:
        raise KeyError("Required column 'Step' not found in CSV")

    # 2) expected components and versions
    versions = ['Simple', 'INC', 'BST']
    # NOTE: 'Insert' is kept here for data cleaning, but excluded from TIME_COMPONENTS for plotting
    components = ['Insert', 'DFS', 'MinI', 'MaxO', 'SBCompute', 'Total']  

    # 3) build canonical column names and check presence (allow small variations)
    # For robustness, accept lowercase variants or missing 'Time' suffix
    def find_col(df_cols, candidates):
        for c in candidates:
            if c in df_cols:
                return c
        return None

    cols = df.columns.tolist()
    for v in versions:
        for comp in components:
            # candidates: e.g., "SimpleInsertTime", "SimpleInsert", "simple_insert_time"
            candidates = [
                f"{v}{comp}Time",
                f"{v}{comp}",
                f"{v}_{comp}Time",
                f"{v}_{comp}",
                f"{v}{comp}time",
                f"{v}_{comp}_time"
            ]
            found = find_col(cols, candidates)
            if not found:
                # Not fatal for some components (you may allow missing SBCompute), but warn/raise as needed
                # Here we'll raise for core components, warn for SuperBubble
                if comp in ['Insert', 'DFS', 'MinI', 'MaxO']:
                    raise KeyError(f"Missing required column for plotting: tried {candidates}")
                else:
                    # Fill missing SBCompute/Total with zeros
                    df[f"{v}{comp}Time"] = 0.0
                    continue
            # rename to canonical form
            canon = f"{v}{comp}Time"
            if found != canon:
                df[canon] = df[found]
    # 4) coerce all canonical time columns to numeric
    for v in versions:
        for comp in components:
            canon = f"{v}{comp}Time"
            if canon in df.columns:
                df[canon] = pd.to_numeric(df[canon], errors='coerce').fillna(0.0)
    return df


def convert_seconds_to_ms(df):
    """Convert all canonical time columns from seconds to milliseconds."""
    df = df.copy()
    versions = ['Simple', 'INC', 'BST']
    components = ['Insert', 'DFS', 'MinI', 'MaxO', 'SBCompute', 'Total']

    for v in versions:
        for comp in components:
            col = f"{v}{comp}Time"
            if col in df.columns:
                df[col] = df[col] * 1000.0
    return df


def smooth_df_rolling(df, window=15):
    """
    Smooth time columns using rolling mean.
    Keeps Step unchanged; only smooths time columns.
    """
    df = df.sort_values("Step").copy()

    versions = ['Simple', 'INC', 'BST']
    components = ['Insert','DFS', 'MinI', 'MaxO', 'SBCompute', 'Total']  # Insert excluded in plots

    for v in versions:
        for comp in components:
            col = f"{v}{comp}Time"
            if col in df.columns:
                df[col] = df[col].rolling(window=window, center=True, min_periods=1).mean()

    return df

def plot_separate_time_components(df, output_dir, dataset_name, scale):
    """
    Generates a 3x2 grid of subplots showing time component performance
    (log scale) vs. Haplotype Step for all three versions.
    Includes a clear legend for the line colors.
    
    --- MODIFIED: Layout changed from 2x3 to 3x2 to accommodate 5 plots ---
    """
    print("Generating Plot 1: Time Component Comparison (Line Plot)...")
    
    # Filter out the 'Base' step which often has zero times and affects log scale visualization
    df_plot = df[df['Step'] != 0].copy()
    
    # Change from (2, 3) to (3, 2) for 5 plots
    fig, axes = plt.subplots(3, 2, figsize=(15, 15))
    axes = axes.flatten() # Flattens the 2D array of axes into a 1D array

    for ax in axes:
        if ax.get_legend() is not None:
            ax.get_legend().remove()

    i = 0
    # Loop through the components (DFS, MinI, MaxO, SBCompute, Total)
    for component, cols in TIME_COMPONENTS.items():
        
        ax = axes[i]
        
        # Plot each version
        eps = 1e-9  # prevents log-scale problems on zeros

        ax.plot(df_plot['Step'], df_plot[cols[0]] + eps, **VERSION_STYLES['Simple'], label='Simple (Red)')
        ax.plot(df_plot['Step'], df_plot[cols[1]] + eps, **VERSION_STYLES['INC'],   label='INC (Blue)')
        ax.plot(df_plot['Step'], df_plot[cols[2]] + eps, **VERSION_STYLES['BST'],   label='BST (Green)')

        ax.set_title(f'{component} Time vs. Haplotype Step — {dataset_name}', fontsize=12)
        if scale == "log":
            ax.set_yscale("log", nonpositive="clip")
        else:
            ax.set_yscale("linear")
        ylabel = f"{component} Time (ms)"
        if scale == "log":
            ylabel += " [Log Scale]"

        ax.set_ylabel(ylabel)
        ax.grid(True, which='both', ls='--', linewidth=0.5, alpha=0.7)
        
        # --- MODIFIED: X-label logic for 3x2 grid ---
        # Only add X-label to bottom row (indices 4 and 5 in the flattened array)
        if i >= 4:
            ax.set_xlabel('Haplotype Step (Ordinal)')
        
        i += 1
        
    # --- ADDED: Turn off the empty 6th subplot ---
    for j in range(i, len(axes)):
        axes[j].axis("off")
        
    # Create the legend based on the plots from the first subplot (axes[0])
    handles, labels = axes[0].get_legend_handles_labels()
    
    legend_handles = [
    Line2D([0], [0], color=VERSION_STYLES['Simple']['color'],
           marker=VERSION_STYLES['Simple']['marker'],
           linestyle=VERSION_STYLES['Simple']['linestyle']),
    Line2D([0], [0], color=VERSION_STYLES['INC']['color'],
           marker=VERSION_STYLES['INC']['marker'],
           linestyle=VERSION_STYLES['INC']['linestyle']),
    Line2D([0], [0], color=VERSION_STYLES['BST']['color'],
           marker=VERSION_STYLES['BST']['marker'],
           linestyle=VERSION_STYLES['BST']['linestyle']),
    ]

    legend_labels = ['Simple (SimpleDFS)', 'INC (IncrementalDFS)', 'BST (BSTIncrementalDFS)']

    # Reserve space at the top and add a single global legend
    fig.subplots_adjust(top=0.92)   # make room for legend
    fig.legend(
        legend_handles,
        legend_labels,
        loc='upper center',
        bbox_to_anchor=(0.5, 0.98),
        ncol=3,
        title='Algorithm (color)'
    )

    # tighten layout (respect the reserved top area) and save including legend
    plt.tight_layout(rect=[0, 0, 1, 0.9])

    # define filename BEFORE saving
    plot_filename = os.path.join(output_dir,f"{dataset_name}_plot1_time_component_comparison.png")
    plt.savefig(plot_filename, dpi=200, bbox_inches='tight')
    print(f"Saved Plot 1 to {plot_filename}")
    plt.close(fig)  # Close the figure


def plot_stacked_bar_breakdown(df, output_dir,dataset_name, num_steps=5):
    """
    Generates a stacked bar chart showing the breakdown of Total Time
    into component times for the Simple, INC, and BST versions across
    a few representative Haplotype Steps.
    """
    print(f"Generating Plot 2: Total Time Stacked Breakdown (First {num_steps} Steps)...")
    
    # Filter out 'Base' and select the first 'num_steps' steps
    df_plot = df[df['Step'] != 0].head(num_steps).copy()
    
    if df_plot.empty:
        print("Not enough non-base data to generate the stacked bar chart.")
        return

    # --- MODIFIED: 'Insert' component has been removed ---
    components = ['Insert','DFS', 'MinI', 'MaxO', 'SBCompute']
    versions = ['Simple', 'INC', 'BST']
    
    # Create the figure
    fig, ax = plt.subplots(figsize=(15, 8))
    
    # Define bar width and positions
    bar_width = 0.25
    step_indices = np.arange(len(df_plot)) # X-coordinates for each step group
    
    # Dictionary to store the cumulative height (bottom) for each version's bar
    version_bottoms = {v: np.zeros(num_steps) for v in versions}
                       
    # X positions for bars: 
    x_positions = {
        'Simple': step_indices - bar_width,
        'INC': step_indices,
        'BST': step_indices + bar_width,
    }

    # Loop through each time component to stack the bars
    for comp in components:
        color = STACK_COLORS[comp]
        
        for version in versions:
            # Column name in the DataFrame, e.g., 'SimpleInsertTime'
            col_name = f'{version}{comp}Time'
            
            # The data for the current component and version across the steps
            data = df_plot[col_name].values

            # Plot the bar segment
            # Only label the component for the 'Simple' version to avoid duplicates in the legend
            label_text = f'{comp} Time' if version == 'Simple' else "" 
            
            ax.bar(
                x_positions[version], 
                data, 
                bar_width, 
                label=label_text,
                bottom=version_bottoms[version], 
                color=color, 
                edgecolor='black',
                linewidth=0.5
            )
            
            # Update the bottom position for the next component in this version's bar
            version_bottoms[version] += data

    # --- Add Version Labels Above Bars and Step Lines ---
    
    # Draw vertical lines to separate step groups
    for x in step_indices[:-1] + (bar_width * 1.5):
        ax.axvline(x=x, color='gray', linestyle=':', linewidth=1.0, alpha=0.5)
    ymax = ax.get_ylim()[1]
    ax.set_ylim(0, ymax * 1.12)
    # Add the version labels (Simple, INC, BST) centered above their respective bars
    # Find the maximum height for each bar group to place the label just above it
    for i, step in enumerate(df_plot['Step']):
        max_heights = {v: version_bottoms[v][i] for v in versions}
        
        for version in versions:
            x_pos = x_positions[version][i]
            y_pos = max_heights[version] + 0.05 * ax.get_ylim()[1] # Offset label slightly above bar
            
            # Place the version label
            ax.text(x_pos, y_pos, version, 
                    ha='center', va='bottom', fontsize=10, 
                    fontweight='bold', color='k')

    # --- Final Plot Aesthetics ---
    
    ax.set_title(f'Total Time Breakdown by Version for First {num_steps} Haplotype Steps — {dataset_name}', fontsize=16)
    ax.set_xlabel('Haplotype Step', fontsize=14)
    ax.set_ylabel('Total Time (ms) [Sum of Components]')
    
    # Set X ticks to be in the center of the step groups
    ax.set_xticks(step_indices) 
    ax.set_xticklabels(df_plot['Step'].astype(str))
    
    # Use a simpler legend focused only on the components (since versions are labeled above bars)
    handles, labels = ax.get_legend_handles_labels()
    # Filter out empty labels from the plotting loop
    valid_handles = [h for h, l in zip(handles, labels) if l]
    valid_labels = [l.replace(' Time', '') for l in labels if l]
    
    ax.legend(valid_handles, valid_labels, loc='upper right', ncol=1, 
              title='Time Component', fancybox=True, shadow=True)
              
    ax.grid(axis='y', ls='--', alpha=0.6)
    plt.tight_layout()
    
    plot_filename = os.path.join(output_dir, f"{dataset_name}_plot2_stacked_time_breakdown.png")
    plt.savefig(plot_filename, dpi=200)
    print(f"Saved Plot 2 to {plot_filename}")
    plt.close(fig)


def plot_haplotype_data(filename, output_dir):
    """
    Reads the haplotype comparison CSV, cleans it using clean_and_normalize_df(),
    and produces three plots:
      1) plot_separate_time_components (3x2 line plots)
      2) plot_aggregated_version_bars (one stacked bar per algorithm)
      3) plot_stacked_bar_breakdown (first N steps stacked bars)

    Usage: plot_haplotype_data("myfile.csv", "plots/")
    """
    # --- Load CSV ---
    try:
        df = pd.read_csv(filename)
        print(f"Successfully loaded '{filename}'. Found {len(df)} rows.")
    except FileNotFoundError:
        print(f"Error: The file '{filename}' was not found.")
        return
    except Exception as e:
        print(f"An error occurred while reading the file '{filename}': {e}")
        return

    # --- Clean & normalize the DataFrame (fails fast with readable errors) ---
    try:
        df = clean_and_normalize_df(df)
    except KeyError as ke:
        print(f"Data cleaning error: {ke}")
        return
    except Exception as e:
        print(f"Unexpected error during data cleaning: {e}")
        return

    # --- Create output directory if it doesn't exist ---
    try:
        os.makedirs(output_dir, exist_ok=True)
    except Exception as e:
        print(f"Error: Could not create output directory '{output_dir}': {e}")
        return

    # --- Plot generation (wrapped so a failure in one plot doesn't stop the others) ---
    try:
        # 1) Component comparison line plots (3x2 subplots)
        plot_separate_time_components(df, output_dir, dataset_name, scale)
    except Exception as e:
        print(f"Error while generating plot_separate_time_components: {e}")

    try:
        # 2) Aggregated single stacked bar per version (Simple/INC/BST)
        plot_aggregated_version_bars(df, output_dir, dataset_name)
    except Exception as e:
        print(f"Error while generating plot_aggregated_version_bars: {e}")

    try:
        # 3) Stacked bar breakdown for the first few steps (default num_steps=5)
        plot_stacked_bar_breakdown(df, output_dir, dataset_name, num_steps=5)
    except Exception as e:
        print(f"Error while generating plot_stacked_bar_breakdown: {e}")

    print(f"\nAll requested plotting attempts finished. Check '{output_dir}' for generated images.")

def plot_aggregated_version_bars(df, output_dir, dataset_name):
    """
    Produce a single stacked bar per version (Simple / INC / BST).
    Each bar is the sum over all steps (excluding Step==0/'Base') of the component times.
    This yields exactly three bars and a clear component legend.
    """
    print("Generating Plot: Aggregated stacked bars (one bar per version)...")
    df_plot = df[df['Step'] != 0].copy()
    if df_plot.empty:
        print("No non-base rows found; skipping aggregated bar plot.")
        return

    versions = ['Simple', 'INC', 'BST']
    # --- MODIFIED: 'Insert' component has been removed ---
    components = ['Insert','DFS', 'MinI', 'MaxO', 'SBCompute']

    # compute sums for each (version, component)
    sums = {v: [] for v in versions}
    for comp in components:
        for v in versions:
            col = f'{v}{comp}Time'
            if col not in df_plot.columns:
                raise KeyError(f"Expected column '{col}' not found.")
            sums[v].append(df_plot[col].sum())

    # plotting
    fig, ax = plt.subplots(figsize=(8, 6))
    x = np.arange(len(versions))
    bottom = np.zeros(len(versions))

    for i, comp in enumerate(components):
        values = [sums[v][i] for v in versions]
        ax.bar(
            x,
            values,
            bottom=bottom,
            width=0.6,
            label=comp,
            color=STACK_COLORS[comp],
            edgecolor='black',
            linewidth=0.4
        )
        bottom += values

    ax.set_xticks(x)
    ax.set_xticklabels(['Simple', 'INC', 'BST'])
    ax.set_ylabel('Time (ms) — summed across steps')
    ax.set_title(f'Aggregated Total Time Breakdown by Version --{dataset_name}')
    ax.legend(title='Component', loc='upper right')
    ax.grid(axis='y', ls='--', alpha=0.5)
    plt.tight_layout()

    plot_filename = os.path.join(output_dir, f"{dataset_name}_plot_aggregated_version_bars.png")
    plt.savefig(plot_filename, dpi=200)
    print(f"Saved aggregated bar plot to {plot_filename}")
    plt.close(fig)


def prepare_peredge_df(df, max_points=5000):
    # Remove Base row
    df = df[df["Step"] != 0].copy()

    # Step already global integer, but enforce
    df["Step"] = pd.to_numeric(df["Step"], errors="coerce").fillna(0).astype(int)
    df = df[df["Step"] > 0].copy()

    # Downsample uniformly if too large (for speed + file size)
    if len(df) > max_points:
        idx = np.linspace(0, len(df) - 1, max_points).astype(int)
        df = df.iloc[idx].copy()

    return df


def plot_separate_time_components_peredge(df, output_dir, dataset_name, scale):
    print("Generating Plot 1 (perEdge): Time Component Comparison (Line Plot)...")

    df_plot = df.copy()

    fig, axes = plt.subplots(3, 2, figsize=(15, 15))
    axes = axes.flatten()

    for ax in axes:
        if ax.get_legend() is not None:
            ax.get_legend().remove()

    i = 0
    for component, cols in TIME_COMPONENTS.items():
        if i >= 6:
            break

        ax = axes[i]
        eps = 1e-9

        ax.plot(df_plot["Step"], df_plot[cols[0]] + eps, **VERSION_STYLES["Simple"])
        ax.plot(df_plot["Step"], df_plot[cols[1]] + eps, **VERSION_STYLES["INC"])
        ax.plot(df_plot["Step"], df_plot[cols[2]] + eps, **VERSION_STYLES["BST"])

        ax.set_title(f"{component} Time vs Edge Updates — {dataset_name}", fontsize=12)
        if scale == "log":
            ax.set_yscale("log", nonpositive="clip")
        else:
            ax.set_yscale("linear")
        ylabel = f"{component} Time (ms)"
        if scale == "log":
            ylabel += " [Log Scale]"  
        ax.set_ylabel(ylabel)
        ax.grid(True, which="both", ls="--", linewidth=0.5, alpha=0.7)

        if i >= 4:
            ax.set_xlabel("Edge Update ID (Global Step)")

        i += 1

    for j in range(i, len(axes)):
        axes[j].axis("off")

    handles = [
        Line2D([0], [0], **VERSION_STYLES["Simple"]),
        Line2D([0], [0], **VERSION_STYLES["INC"]),
        Line2D([0], [0], **VERSION_STYLES["BST"]),
    ]
    labels = ["Simple", "INC", "BST"]

    fig.subplots_adjust(top=0.92)
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.98), ncol=3)

    plt.tight_layout(rect=[0, 0, 1, 0.9])
    plot_filename = os.path.join(output_dir, f"{dataset_name}_plot1_time_component_comparison_peredge.png")
    plt.savefig(plot_filename, dpi=200, bbox_inches="tight")
    print(f"Saved Plot 1 to {plot_filename}")
    plt.close(fig)


def plot_stacked_bar_breakdown_peredge(df, output_dir, dataset_name, num_steps=10):
    print(f"Generating Plot 2 (perEdge): Total Time Stacked Breakdown (First {num_steps} Edge Updates)...")

    df_plot = df.head(num_steps).copy()
    if df_plot.empty:
        print("Not enough data to generate stacked bar chart.")
        return

    components = ['Insert',"DFS", "MinI", "MaxO", "SBCompute"]
    versions = ["Simple", "INC", "BST"]

    fig, ax = plt.subplots(figsize=(15, 8))

    bar_width = 0.25
    step_indices = np.arange(len(df_plot))

    version_bottoms = {v: np.zeros(len(df_plot)) for v in versions}

    x_positions = {
        "Simple": step_indices - bar_width,
        "INC": step_indices,
        "BST": step_indices + bar_width,
    }

    # stack bars
    for comp in components:
        color = STACK_COLORS[comp]
        for version in versions:
            col_name = f"{version}{comp}Time"
            data_vals = df_plot[col_name].values

            label_text = f"{comp}" if version == "Simple" else ""

            ax.bar(
                x_positions[version],
                data_vals,
                bar_width,
                label=label_text,
                bottom=version_bottoms[version],
                color=color,
                edgecolor="black",
                linewidth=0.5,
            )

            version_bottoms[version] += data_vals

    # vertical separators between edge groups (like perbatch plot)
    for x in step_indices[:-1] + (bar_width * 1.5):
        ax.axvline(x=x, color="gray", linestyle=":", linewidth=1.0, alpha=0.5)

    # ✅ add Simple/INC/BST labels on top of each bar
    for i in range(len(df_plot)):
        for version in versions:
            x_pos = x_positions[version][i]
            y_pos = version_bottoms[version][i] + 0.03 * ax.get_ylim()[1]  # small offset
            ax.text(
                x_pos, y_pos, version,
                ha="center", va="bottom",
                fontsize=10, fontweight="bold", color="black"
            )

    ax.set_xticks(step_indices)
    ax.set_xticklabels(df_plot["Step"].astype(str), rotation=0)

    ax.set_xlabel("Edge Update ID (Global Step)")
    ax.set_ylabel("Time (ms)")
    ax.set_title(f"Total Time Breakdown by Version for First {num_steps} Edge Updates -- {dataset_name}")

    ax.legend(title="Component", loc="upper right")
    ax.grid(True, axis="y", linestyle="--", alpha=0.5)

    plot_filename = os.path.join(output_dir, f"{dataset_name}_plot2_stacked_total_breakdown_peredge.png")
    plt.tight_layout()
    plt.savefig(plot_filename, dpi=200, bbox_inches="tight")
    print(f"Saved Plot 2 to {plot_filename}")
    plt.close(fig)

def run_perbatch(csv_path, output_dir,dataset_name, scale):
    df = pd.read_csv(csv_path)
    df = clean_and_normalize_df(df)
    # convert to ms
    df = convert_seconds_to_ms(df)

    # smoothing (batch steps are fewer, so small window)
    df_smooth = smooth_df_rolling(df, window=3)
    plot_separate_time_components(df_smooth, output_dir,dataset_name, scale)
    plot_stacked_bar_breakdown(df, output_dir,dataset_name, num_steps=5)
    plot_aggregated_version_bars(df, output_dir,dataset_name)
    
    plot_total_time_bar(df, output_dir, dataset_name)


def run_peredge(csv_path, output_dir,dataset_name, scale):
    df = pd.read_csv(csv_path)
    df = clean_and_normalize_df(df)

    # remove Base and sort
    df = df[df["Step"] != 0].copy()
    df = df.sort_values("Step")

    # convert to ms
    df = convert_seconds_to_ms(df)

    # OPTIONAL downsample (still fine)
    max_points = 5000
    if len(df) > max_points:
        idx = np.linspace(0, len(df) - 1, max_points).astype(int)
        df = df.iloc[idx].copy()

    # smoothing for edge updates (bigger window)
    df_smooth = smooth_df_rolling(df, window=50)

    plot_separate_time_components_peredge(df_smooth, output_dir,dataset_name, scale)
    plot_stacked_bar_breakdown_peredge(df, output_dir,dataset_name, num_steps=10)
    plot_aggregated_version_bars(df, output_dir,dataset_name)
    
    plot_total_time_bar(df, output_dir, dataset_name)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True, help="Path to CSV file generated by parse_to_csv.py")
    parser.add_argument("--out", required=True, help="Output directory for plots")
    parser.add_argument("--mode", choices=["perbatch", "peredge"], required=True)
    parser.add_argument("--dataset", required=True, help="Dataset name")
    parser.add_argument(
    "--scale",
    choices=["linear", "log"],
    default="linear",
    help="Y-axis scale for line plots"
    )


    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    if args.mode == "perbatch":
        run_perbatch(args.csv, args.out, args.dataset, args.scale)
    else:
        run_peredge(args.csv, args.out,  args.dataset, args.scale)

    print("✅ Plot generation complete!")

if __name__ == "__main__":

    main()
    """
    Main entry point for the script.
    Expects two command-line arguments:
    1. The path to the input CSV file.
    2. The path to the output directory for plots.
    
    Example:
    python3 plotcsv.py --csv file.csv --out plots --mode perbatch --dataset Ecoli
    python3 plotcsv.py --csv file.csv --out plots --mode peredge --dataset Ecoli
    python3 plotcsv.py --csv file.csv --out plots --mode peredge --dataset Ecoli --scale linear
    python3 plotcsv.py --csv file.csv --out plots --mode peredge --dataset Ecoli --scale log

    """
    
    # Ensure the user has the required packages
    # You can install them with:
    # pip install pandas matplotlib numpy
    """
    # Check if a filename AND output dir were provided as command-line arguments
    if len(sys.argv) > 2:
        # sys.argv[0] is the script name, sys.argv[1] is the input file, sys.argv[2] is the output dir
        filepath = sys.argv[1]
        output_path = sys.argv[2]
        plot_haplotype_data(filepath, output_path)
    else:
        # If not enough arguments are provided, show the user how to run the script
        print("Usage: python plotcsv.py <path_to_your_csv_file> <path_to_output_directory>")
        print(f"Example: python plotcsv.py haplotype_build_comparison.csv plots/")
    """    