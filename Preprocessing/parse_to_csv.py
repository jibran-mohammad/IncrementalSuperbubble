import pandas as pd
import re
import os
import pathlib
import sys
import argparse
# ==============================================================================
#  Dynamic Path Discovery Function (from graphRunner.py)
# ==============================================================================

def find_repo_root(start_path: pathlib.Path, repo_name: str) -> pathlib.Path:
    """
    Finds the absolute path to a named repository root directory.
    
    It starts from 'start_path' (usually the script's directory)
    and walks up the filesystem tree until it finds a directory
    matching 'repo_name'.
    """
    # Start from the resolved, absolute path
    current_dir = start_path.resolve()

    # Loop until we hit the filesystem root (where parent == self)
    while current_dir != current_dir.parent:
        if current_dir.name == repo_name:
            return current_dir  # Found it!
        current_dir = current_dir.parent

    # Final check in case the root directory itself is the repo_name
    if current_dir.name == repo_name:
        return current_dir

    # If we get here, we've hit the filesystem root and found nothing
    raise FileNotFoundError(
        f"Could not find repository root '{repo_name}' in any parent path of {start_path}"
    )

# ==============================================================================
#  Helper Function for Time Conversion
# ==============================================================================

def robust_float_convert(val_str: str) -> float:
    """
    Converts a string to a float, correctly handling
    standard decimal or scientific notation (e.g., "1.23e-05").
    Returns 0.0 on conversion error.
    """
    try:
        return float(val_str)
    except ValueError:
        print(f"    Warning: Could not convert '{val_str}' to float. Using 0.0.")
        return 0.0
    except TypeError:
        print(f"    Warning: Invalid type for float conversion. Using 0.0.")
        return 0.0

# ==============================================================================
#  Parsing Function (Updated)
# ==============================================================================

def parse_file(filename, time_column_name):
    """
    Parses one of the superbubble output files.
    
    The file has a unique format where graph stats (edges, nodes)
    are on one line, and the build stats (name, time, superbubbles)
    are on the next line.
    """
    data = []
    
    # Check if file exists before trying to open
    # This now works with pathlib.Path objects
    if not os.path.exists(filename):
        print(f"    Warning: File not found: {os.path.basename(filename)}. Skipping.")
        return pd.DataFrame()

    with open(filename, 'r') as f:
        lines = f.readlines()

    # Use an iterator to process lines in pairs
    line_iter = iter(lines)
    
    try:
        while True:
            # --- Line 1: Edges and Nodes ---
            # Example: "8249 4203"
            stats_line = next(line_iter).strip()
            # Stop if we hit the peak memory line
            if "PeakMemory" in stats_line:
                break
                
            stats_parts = stats_line.split()
            if len(stats_parts) < 2:
                continue # Skip malformed lines
                
            edges = int(stats_parts[0])
            nodes = int(stats_parts[1])

            # --- Line 2: Name, Time, Superbubbles ---
            # Example 1: "BaseGraph_Build 0.00223163 3926"
            # Example 2: "1 GCF_000618905.3... 0.0311279 4192"
            build_line = next(line_iter).strip()
            build_parts = build_line.split()
            
            if not build_parts:
                continue # Skip empty lines

            step_name = ""
            time_val = 0.0
            superbubbles = 0
            
            if build_parts[0] == "BaseGraph_Build":
                step = 0  # Use 0 for BaseGraph for sorting
                step_name = "BaseGraph_Build"
                # Use robust conversion
                time_val = robust_float_convert(build_parts[1])
                superbubbles = int(build_parts[2])
            else:
                # Handle lines with a step number
                step = int(build_parts[0]) # Convert step to int
                step_name = build_parts[1]
                # Use robust conversion
                time_val = robust_float_convert(build_parts[2])
                superbubbles = int(build_parts[3])

            data.append({
                'Step': step,
                'Haplotype_Name': step_name,
                'Nodes': nodes,
                'Edges': edges,
                'Superbubbles': superbubbles,
                time_column_name: time_val
            })
            
    except StopIteration:
        # End of file
        pass
    except Exception as e:
        print(f"    Error parsing {filename}: {e}")

    return pd.DataFrame(data)

# ==============================================================================
#  New Processing Function
# ==============================================================================

def process_time_directory(time_dir_path: pathlib.Path):
    """
    Parses all three time files within a specific 'timeFile' directory
    and saves the combined CSV inside that same directory.
    """
    print(f"\nProcessing directory: {time_dir_path}")
    
    if not time_dir_path.is_dir():
        print(f"  Warning: Directory not found. Skipping.")
        return

    # Define the files and their corresponding time column names
    # We now use the full path by joining the directory
    files_to_parse = {
        time_dir_path / "SuperbubbleSimpleTimeOutput.txt": "Simple_Time_s",
        time_dir_path / "SuperbubbleIncrementalDFSOutput.txt": "DFS_Time_s",
        time_dir_path / "SuperbubbleIncrementalAVLDFSOutput.txt": "AVLDFS_Time_s"
    }

    all_data = []
    
    # Parse each file
    for filepath, time_col in files_to_parse.items():
        df = parse_file(filepath, time_col)
        if not df.empty:
            all_data.append(df)

    if not all_data:
        print("  No data was parsed. No CSV will beG be generated.")
        return

    # --- Merge the data ---
    try:
        merged_df = all_data[0].copy()
        merged_df = merged_df.set_index('Step')

        for df in all_data[1:]:
            time_col_name = [col for col in df.columns if col.endswith('_Time_s')][0]
            other_df_time = df.set_index('Step')[[time_col_name]]
            merged_df = merged_df.join(other_df_time, how='outer')

        merged_df = merged_df.reset_index()
        
        # --- NEW: Sort by the integer 'Step' column ---
        merged_df = merged_df.sort_values(by='Step').reset_index(drop=True)

        # --- NEW: Convert 'Step' column back to desired string format ---
        # Change type to 'object' to allow storing a mix of strings ("Base") and numbers
        merged_df['Step'] = merged_df['Step'].astype(object)
        # Find the row where Step is 0 (our placeholder for Base) and replace it
        merged_df.loc[merged_df['Step'] == 0, 'Step'] = 'Base'

        common_cols = ['Step', 'Haplotype_Name', 'Nodes', 'Edges', 'Superbubbles']
        time_cols = [col for col in merged_df.columns if col.endswith('_Time_s')]
        final_columns = common_cols + sorted(time_cols)
        
        final_columns = [col for col in final_columns if col in merged_df.columns]
        merged_df = merged_df[final_columns]

        # --- Save to CSV ---
        output_filename = time_dir_path / "haplotype_build_comparison.csv"
        merged_df.to_csv(output_filename, index=False)
        
        print(f"  ✅ Successfully saved data to: {output_filename}")

    except Exception as e:
        print(f"  ❌ Error merging or saving data for {time_dir_path}: {e}")


# ==============================================================================
#  New Dynamic Main Function
# ==============================================================================

def main():
    """
    Main function to find output directories and process them.

    --------------------------------------------------------------------------
    HOW TO RUN THIS SCRIPT:
    --------------------------------------------------------------------------

    1.  **Prerequisites:**
        * You must have the 'pandas' library installed:
            `pip install pandas`
        * This script MUST be located somewhere inside your
            'IncSuperbubble' repository directory structure.
            It doesn't matter where (e.g., can be in 'Preprocessing',
            'output', or the root), as it will find the root dynamically.

    2.  **Running the script:**
        * Open your terminal or command prompt.
        * Navigate to the directory containing this script.
        * Run the script using python, specifying the dataset:

            `python3 parse_to_csv.py -d ecoli`
            (or `python3 parse_to_csv.py --dataset human`)
            (or `python3 parse_to_csv.py -d all`)

    3.  **What it does:**
        * It automatically finds the 'IncSuperbubble' repository root.
        * It looks for the 'output' directory.
        * It processes `output/Ecoli/timeFile/` and saves a CSV inside it.
        * It finds all `output/Human/chromosome_*` directories, looks
            for a `timeFile` subdirectory in each one, and saves a
            separate CSV inside each `timeFile` directory it finds.

    4.  **Output:**
        * A file named `haplotype_build_comparison.csv` will be created
            inside each `timeFile` directory that is processed.
    """
    parser = argparse.ArgumentParser(
        description="Parses haplotype time-tracking files into CSVs.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-d", "--dataset",
        choices=["ecoli", "human", "all"],
        required=True,
        help="The dataset to process. 'ecoli', 'human', or 'all'."
    )
    args = parser.parse_args()

    try:
        # 1. Define the name of the repository we are looking for
        REPO_NAME = "IncSuperbubble"
        
        # 2. Get the directory containing this script
        try:
            # Standard way (when running as script)
            script_dir = pathlib.Path(__file__).resolve().parent
        except NameError:
            # Fallback (for interactive, etc.)
            script_dir = pathlib.Path.cwd()

        # 3. Find the repo root path object
        REPO_ROOT_PATH = find_repo_root(script_dir, REPO_NAME)
        OUTPUT_DIR = REPO_ROOT_PATH / "output"
        
        print(f"ℹ️ Located repository root: {REPO_ROOT_PATH}")
        print(f"ℹ️ Looking for output in: {OUTPUT_DIR}")

    except FileNotFoundError as e:
        print(f"❌ Critical Error: {e}", file=sys.stderr)
        print("Please ensure this script is located somewhere inside the 'IncSuperbubble' repository.", file=sys.stderr)
        sys.exit(1)

    # --- 1. Process Ecoli Dataset ---
    if args.dataset in ["ecoli", "all"]:
        ecoli_time_dir = OUTPUT_DIR / "Ecoli" / "timeFile"
        process_time_directory(ecoli_time_dir)

    # --- 2. Process Human Dataset ---
    if args.dataset in ["human", "all"]:
        human_output_dir = OUTPUT_DIR / "Human"
        if not human_output_dir.is_dir():
            print(f"\nℹ️ 'Human' directory not found at {human_output_dir}, skipping.")
        else:
            print(f"\nScanning {human_output_dir} for chromosome directories...")
            # Find all directories starting with "chromosome_"
            chromosome_dirs = sorted([d for d in human_output_dir.iterdir() if d.is_dir() and d.name.startswith("chromosome_")])
            
            if not chromosome_dirs:
                print("  No 'chromosome_*' directories found.")
                # This 'pass' is fine, or we can just let it continue to the end
            else:
                print(f"  Found {len(chromosome_dirs)} chromosome directories to process.")
                
                # Loop over each chromosome directory and process its 'timeFile' subdir
                for chrom_dir in chromosome_dirs:
                    human_time_dir = chrom_dir / "timeFile"
                    process_time_directory(human_time_dir)

    print("\n🎉 All processing complete!")


if __name__ == "__main__":
    main()