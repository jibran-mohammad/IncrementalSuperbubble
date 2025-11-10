#!/usr/bin/env python3

"""
Parses all 'SuperbubbleSimpleTimeOutput.txt' style files in a specified
input folder and generates a corresponding .tex file for each in an
output folder.
"""

import sys
import os # <-- Added for path handling
from pathlib import Path # <-- Added for modern path handling

def get_latex_preamble():
# ... existing code ...
    return r"""
\begin{document}

\maketitle
\section*{Summary}
This report details the performance metrics for the superbubble computation process.
"""

def get_latex_document_end():
# ... existing code ...
def parse_data(input_file):
    """Parses the input text file and returns structured data."""
    try:
        with open(input_file, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.", file=sys.stderr)
        return None, None, None # <-- Changed from sys.exit
    except Exception as e:
        print(f"Error reading file '{input_file}': {e}", file=sys.stderr)
        return None, None, None # <-- Changed from sys.exit

    if len(lines) < 5:
        print(f"Error: File '{input_file}' appears to be too short or in an unexpected format.", file=sys.stderr)
        return None, None, None # <-- Changed from sys.exit

    # --- Parse Base Graph (Lines 0, 1) ---
    try:
# ... existing code ...
        base_graph_info = {
            "nodes": base_nodes,
            "edges": base_edges,
            "build_time": base_build_time
        }
    except Exception as e:
        print(f"Error parsing base graph in '{input_file}': {e}", file=sys.stderr)
        return None, None, None

    # --- Parse Peak Memory (Last Line) ---
    peak_memory_line = lines[-1]
# ... existing code ...
    if "PeakMemory(KB):" in peak_memory_line:
        peak_memory = peak_memory_line.split()[-1]

    # --- Parse Haplotype Blocks (Middle Lines) ---
# ... existing code ...
    return base_graph_info, haplotypes, peak_memory


def generate_latex_output(base_info, haplotypes, peak_memory, input_file_name):
    """Generates the full LaTeX string from the parsed data."""
    
    output = []
    # --- Pass the original filename to the \author command ---
    preamble = get_latex_preamble().replace(
        r'\author{Data from \texttt{SuperbubbleSimpleTimeOutput.txt}}',
        fr'\author{{Data from \texttt{{{input_file_name.replace("_", r"\_")}}}}}'
    )
    output.append(preamble)
    output.append(get_latex_document_start())

    # --- Add Base Graph Section ---
# ... existing code ...
    return "\n".join(output)

def main():
    """Main function to parse arguments and run the script."""
    
    # --- Define FOLDER paths here ---
    # !! Replace these with the full, absolute paths on your system !!
    input_folder_path = Path("/home/dell/Research/Superbubble/Experiment/IncSuperbubble/output/Ecoli/timeFile")
    output_folder_path = Path("/home/dell/Research/Superbubble/Experiment/IncSuperbubble/output/Ecoli/timeFile/table")
    
    # --- Create output folder if it doesn't exist ---
    try:
        output_folder_path.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        print(f"Error creating output directory '{output_folder_path}': {e}", file=sys.stderr)
        sys.exit(1)

    # --- Find and process all .txt files in the input folder ---
    txt_files = list(input_folder_path.glob('*.txt'))
    
    if not txt_files:
        print(f"No .txt files found in '{input_folder_path}'.")
        return

    print(f"Found {len(txt_files)} .txt file(s) to process.")
    print(f"Output will be saved to '{output_folder_path}'")
    
    processed_count = 0
    error_count = 0

    for input_file in txt_files:
        print(f"\n--- Processing: {input_file.name} ---")
        
        base_info, haplotypes, peak_memory = parse_data(input_file)
        
        # --- Skip file if parsing failed ---
        if base_info is None:
            print(f"Skipping file '{input_file.name}' due to parsing error.")
            error_count += 1
            continue

        print(f"  - Found Base Graph info.")
        print(f"  - Found {len(haplotypes)} haplotype entries.")
        print(f"  - Found Peak Memory: {peak_memory} KB.")
        
        print(f"  - Generating LaTeX code...")
        # Pass the input file name to customize the title
        latex_content = generate_latex_output(base_info, haplotypes, peak_memory, input_file.name)
        
        # --- Define the output .tex file path ---
        output_file = output_folder_path / input_file.with_suffix('.tex').name
        
        try:
            with open(output_file, 'w') as f:
                f.write(latex_content)
            print(f"  - Success! LaTeX file saved to '{output_file}'")
            processed_count += 1
        except Exception as e:
            print(f"  - Error writing to output file '{output_file}': {e}", file=sys.stderr)
            error_count += 1

    print(f"\n--- Batch processing complete ---")
    print(f"Successfully processed: {processed_count} file(s)")
    print(f"Failed/