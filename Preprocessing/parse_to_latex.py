#!/usr/bin/env python3

"""
Parses multiple algorithm output files from an input folder.
Generates:
1. A single, comprehensive LaTeX (.tex) report for COMPARISON.
2. Separate detailed .tex files for EACH algorithm's full table.
"""

import sys
import re
import datetime  # --- FIXED: Added missing import ---
from pathlib import Path
from decimal import Decimal # For precise handling of scientific notation

# --- CONFIGURATION ---

# **IMPORTANT**: Update this path to your input folder
INPUT_FOLDER = Path("/home/dell/Research/Superbubble/Experiment/IncSuperbubble/output/Ecoli/timeFile")

# **NEW**: Specify a separate folder for your final .tex report
OUTPUT_FOLDER = Path("/home/dell/Research/Superbubble/Experiment/IncSuperbubble/output/Ecoli/timeFile/table")

# **IMPORTANT**: This is the mapping of file prefixes to "pretty names"
ALGO_NAME_MAP = {
    # Example mapping:
    "SuperbubbleSimpleTimeOutput": "Static",
    "SuperbubbleIncrementalDFSOutput": "Dynamic",
    "SuperbubbleIncrementalAVLDFSOutput": "DynamicBST",
    # Add your actual algorithm file prefixes and names here
    # "ecoli_inc_output": "Incremental Algorithm",
    # "ecoli_static_output": "Static Algorithm",
}

# **NEW**: The name of the main comparison report file
COMPARISON_REPORT_FILE = "ecoli_algorithm_comparison.tex"

# --- END CONFIGURATION ---


def parse_time(time_str: str) -> float:
    """
    Parses a time string, handling scientific notation.
    Returns time in seconds.
    """
    try:
        # Use Decimal for precision, then convert to float
        return float(Decimal(time_str))
    except Exception:
        print(f"  [Warning] Could not parse time '{time_str}', using 0.0", file=sys.stderr)
        return 0.0


# --- MODIFIED: This entire function is rewritten to correctly parse the 2-line block format ---
def parse_data_file(file_path: Path) -> dict:
    """
    Parses a single algorithm output file.
    
    File Format Expected (2-line blocks):
    Line 1: <edges> <nodes>
    Line 2: <(id) name> <time> <superbubbles>
    ... (blocks repeat) ...
    
    Also handles split <edges> <nodes> lines.
    Last Line: PeakMemory(KB): <memory>
    """
    try:
        with open(file_path, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]
    except Exception as e:
        print(f"  [Error] Could not read file '{file_path}': {e}", file=sys.stderr)
        return {}

    parsed_data = {
        "base_graph": None,
        "haplotypes": [],
        "peak_memory_kb": None
    }

    i = 0
    while i < len(lines):
        line = lines[i]

        # --- Check 1: Is this the Peak Memory line? ---
        memory_match = re.search(r"PeakMemory\(KB\):\s*(\d+)", line, re.IGNORECASE)
        if memory_match:
            parsed_data["peak_memory_kb"] = int(memory_match.group(1))
            i += 1
            continue

        # --- Check 2: Is this a data block? ---
        # This line should be <edges> <nodes>
        edges = 0
        nodes = 0
        data_line_index = i + 1
        increment = 2 # Default for 2-line block

        edges_nodes_match = re.match(r"^\s*(\d+)\s+(\d+)\s*$", line)
        
        if edges_nodes_match:
            edges = int(edges_nodes_match.group(1))
            nodes = int(edges_nodes_match.group(2))
        elif i + 1 < len(lines):
            # Check for split-line case (e.g., "10952" and "6906")
            line1_match = re.match(r"^\s*(\d+)\s*$", lines[i])
            line2_match = re.match(r"^\s*(\d+)\s*$", lines[i+1])
            if line1_match and line2_match:
                edges = int(line1_match.group(1))
                nodes = int(line2_match.group(1))
                data_line_index = i + 2
                increment = 3 # Consumed 3 lines (edges, nodes, data)
            else:
                # Unrecognized line
                i += 1
                continue
        else:
            # Unrecognized line / single line at end of file
            i += 1
            continue
            
        # Check if data line exists
        if data_line_index >= len(lines):
            i += 1 # End of file
            continue
            
        # --- Parse the Data Line ---
        try:
            data_line = lines[data_line_index]
            parts = data_line.split()
            
            if len(parts) < 3: # Must have at least name, time, superbubbles
                i += increment # Skip malformed block
                continue
                
            superbubbles = int(parts[-1])
            time_val = parse_time(parts[-2])
            identifier = parts[0]
            
            block_data = {
                "nodes": nodes,
                "edges": edges,
                "time": time_val,
                "superbubbles": superbubbles
            }

            if identifier.isdigit():
                # It's a Haplotype: "1 GCF_... 0.01 4192"
                block_data["id"] = int(identifier)
                block_data["name"] = " ".join(parts[1:-2]) # Name is everything between id and time
                parsed_data["haplotypes"].append(block_data)
            else:
                # It's BaseGraph: "BaseGraph_Build 0.002 3926"
                block_data["name"] = " ".join(parts[0:-2]) # Name is everything before time
                if "basegraph" in block_data["name"].lower():
                    parsed_data["base_graph"] = block_data
            
            i += increment # Move to the next block

        except Exception as e:
            # Failed to parse this block, move to the next line to be safe
            # print(f"  [Debug] Failed to parse block at line {i+1}: {e}", file=sys.stderr)
            i += 1

    # --- Post-processing and Sanity Checks ---
    if not parsed_data["base_graph"] and not parsed_data["haplotypes"]:
        print(f"  [Error] No valid BaseGraph or Haplotype blocks found in '{file_path.name}'.", file=sys.stderr)
        return {}

    # Calculate totals
    total_time = 0.0
    total_haplo_time = 0.0
    
    if parsed_data["base_graph"]:
        total_time += parsed_data["base_graph"].get("time", 0.0)

    for hap in parsed_data["haplotypes"]:
        hap_time = hap.get("time", 0.0)
        total_time += hap_time
        total_haplo_time += hap_time

    parsed_data["total_time"] = total_time
    parsed_data["total_haplotype_time"] = total_haplo_time
    
    if parsed_data["haplotypes"]:
        parsed_data["avg_haplotype_time"] = total_haplo_time / len(parsed_data["haplotypes"])
    else:
        parsed_data["avg_haplotype_time"] = 0.0

    return parsed_data
# --- END OF MODIFIED FUNCTION ---


def make_tex_safe(text: str) -> str:
    """Replaces characters that have special meaning in LaTeX."""
    return text.replace("_", r"\_").replace("%", r"\%").replace("&", r"\&")


def generate_comparison_fragment(all_data: dict) -> str:
    """
    Generates the LaTeX code for the main comparison summary table
    and the performance analysis text.
    
    Returns a LaTeX string (fragment).
    """
    output = []
    output.append(r"\section*{Comparison Summary}")
    output.append(r"% --- Generated by script on " + f"{datetime.datetime.now()} ---")
    
    # --- 1. Comparison Table ---
    output.append(r"\begin{table}[ht]")
    output.append(r"\centering")
    output.append(r"\caption{Algorithm Performance Summary}")
    output.append(r"\label{tab:comparison_summary}")
    output.append(r"\begin{tabular}{l r r r r}")
    output.append(r"\toprule") # --- REVERTED TO \toprule ---
    output.append(r"\textbf{Algorithm} & \textbf{Total Time (s)} & \textbf{BaseGraph Time (s)} & \textbf{Total Haplo. Time (s)} & \textbf{Avg. Haplo. Time (s)} \\")
    output.append(r"\midrule") # --- REVERTED TO \midrule ---

    # Sort data by total time, fastest first
    sorted_algos = sorted(all_data.items(), key=lambda item: item[1].get("total_time", float('inf')))

    for algo_name, data in sorted_algos:
        safe_name = make_tex_safe(algo_name)
        base_time = data.get("base_graph", {}).get("time", 0.0)
        output.append(
            f"{safe_name} & "
            f"{data.get('total_time', 0.0):.6f} & "
            f"{base_time:.6f} & "
            f"{data.get('total_haplotype_time', 0.0):.6f} & "
            f"{data.get('avg_haplotype_time', 0.0):.6f} \\\\"
        )
    
    output.append(r"\bottomrule") # --- REVERTED TO \bottomrule ---
    output.append(r"\end{tabular}")
    output.append(r"\end{table}")
    
    # --- 2. Performance Analysis Text ---
    output.append(r"\subsection*{Performance Analysis}")
    if len(sorted_algos) < 2:
        output.append("Only one algorithm was processed, so no comparison is available.")
        return "\n".join(output)

    # Compare everything to the fastest algorithm
    fastest_name = make_tex_safe(sorted_algos[0][0])
    fastest_time = sorted_algos[0][1].get("total_time", 0.0)

    comparisons = []
    if fastest_time > 0:
        for algo_name, data in sorted_algos[1:]:
            current_time = data.get("total_time", 0.0)
            if current_time > 0:
                speed_factor = current_time / fastest_time
                comparisons.append(
                    f"Overall, \\textbf{{{fastest_name}}} was \\textbf{{{speed_factor:.2f}x faster}} than \\textbf{{{make_tex_safe(algo_name)}}}"
                    f" ({fastest_time:.4f}s vs {current_time:.4f}s)."
                )

    if comparisons:
        output.append(r"\begin{itemize}")
        for comp in comparisons:
            output.append(rf"    \item {comp}")
        output.append(r"\end{itemize}")
    else:
        output.append("All algorithms had similar performance or insufficient data for comparison.")

    return "\n".join(output)


def generate_algo_table_fragment(algo_name: str, data: dict) -> str:
    """
    Generates the LaTeX code for a single algorithm's detailed
    longtable, including BaseGraph and all haplotypes.
    
    Returns a LaTeX string (fragment).
    """
    safe_name = make_tex_safe(algo_name)
    output = []
    output.append(r"\subsection*{" + f"Detailed Metrics: {safe_name}" + r"}")
    output.append(r"% --- Generated by script on " + f"{datetime.datetime.now()} ---")
    
    # Use \begin{small} to make the table fit better
    output.append(r"\begin{small}")
    output.append(r"\begin{longtable}{l r r r r}")

    # ---
    # THIS IS THE CORRECT HEADER STRUCTURE FOR longtable
    # ---
    # 1. Caption
    output.append(r"\caption{Detailed metrics for " + safe_name + r"} \label{tab:" + re.sub(r'[^a-zA-Z0-9]', '', algo_name) + r"} \\")
    
    # 2. Header for first page
    output.append(r"\toprule") # --- REVERTED TO \toprule ---
    output.append(r"\textbf{Component} & \textbf{Nodes} & \textbf{Edges} & \textbf{Time (s)} & \textbf{Superbubbles} \\")
    output.append(r"\midrule") # --- REVERTED TO \midrule ---
    output.append(r"\endfirsthead")
    
    # 3. Header for all continued pages
    output.append(r"\caption[]{ (continued)} \\")
    output.append(r"\toprule") # --- REVERTED TO \toprule ---
    output.append(r"\textbf{Component} & \textbf{Nodes} & \textbf{Edges} & \textbf{Time (s)} & \textbf{Superbubbles} \\")
    output.append(r"\midrule") # --- REVERTED TO \midrule ---
    output.append(r"\endhead")
    
    # 4. Footer for all pages except the last
    output.append(r"\bottomrule") # --- REVERTED TO \bottomrule ---
    output.append(r"\multicolumn{5}{r}{\textit{Continued on next page...}} \\")
    output.append(r"\endfoot")

    # 5. Footer for the very last page
    output.append(r"\bottomrule") # --- REVERTED TO \bottomrule ---
    output.append(r"\endlastfoot")
    # ---
    # END OF HEADER/FOOTER DEFINITIONS
    # ---

    # --- Table Body ---
    
    # Add Base Graph
    if "base_graph" in data and data["base_graph"]:
        bg = data["base_graph"]
        output.append(
            rf"\textbf{{{make_tex_safe(bg['name'])}}} & "
            f"{bg['nodes']} & {bg['edges']} & {bg['time']:.6f} & {bg['superbubbles']} \\\\"
        )
    else:
        output.append(r"\textbf{BaseGraph} & N/A & N/A & N/A & N/A \\")

    # Add Haplotypes
    for hap in data.get("haplotypes", []):
        output.append(
            f"Haplotype {hap['id']} & "
            f"{hap['nodes']} & {hap['edges']} & {hap['time']:.6f} & {hap['superbubbles']} \\\\"
        )

    # --- Add Totals Row ---
    output.append(r"\midrule") # --- REVERTED TO \midrule ---
    output.append(
        rf"\textbf{{Total (Base + Haplotypes)}} & & & "
        rf"\textbf{{{data.get('total_time', 0.0):.6f}}} & \\"
    )
    
    if "peak_memory_kb" in data and data["peak_memory_kb"]:
         output.append(
            rf"\textbf{{Peak Memory (KB)}} & & & & \textbf{{{data['peak_memory_kb']}}} \\"
        )
    
    # --- End Table ---
    output.append(r"\end{longtable}")
    output.append(r"\end{small}")
    
    return "\n".join(output)


def main():
    """Main execution function."""
    print("Starting LaTeX Report Generation...")
    print(f"Input Folder: {INPUT_FOLDER}")
    print(f"Output Folder: {OUTPUT_FOLDER}")

    # --- Create output folder ---
    try:
        OUTPUT_FOLDER.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        print(f"\n[ERROR] Could not create output folder '{OUTPUT_FOLDER}': {e}", file=sys.stderr)
        sys.exit(1)

    all_algo_data = {}
    print("\n--- Parsing Algorithm Outputs ---")

    for file_prefix, algo_name in ALGO_NAME_MAP.items():
        print(f"Looking for algorithm: '{algo_name}' (prefix: '{file_prefix}*')")
        
        # Find the first file that matches the prefix
        found_files = list(INPUT_FOLDER.glob(f"{file_prefix}*"))
        if not found_files:
            print(f"  [Warning] No file found with prefix '{file_prefix}' in '{INPUT_FOLDER}'", file=sys.stderr)
            continue
            
        file_to_parse = found_files[0]
        if len(found_files) > 1:
            print(f"  [Warning] Found {len(found_files)} matches; using the first one: '{file_to_parse.name}'", file=sys.stderr)

        print(f"  - Parsing {file_to_parse.name}...")
        parsed_data = parse_data_file(file_to_parse)

        if parsed_data:
            all_algo_data[algo_name] = parsed_data
            print(f"  - Success: Parsed BaseGraph and {len(parsed_data['haplotypes'])} haplotypes.")
        else:
            print(f"  [Error] Failed to parse data for '{algo_name}' from {file_to_parse.name}", file=sys.stderr)

    if not all_algo_data:
        print("\n[ERROR] No valid data was parsed from any file. Aborting.", file=sys.stderr)
        sys.exit(1)

    print("\n--- Generating LaTeX Fragments ---")

    # --- 1. Generate and save the main comparison file ---
    try:
        comparison_latex = generate_comparison_fragment(all_algo_data)
        comparison_file_path = OUTPUT_FOLDER / COMPARISON_REPORT_FILE
        
        with open(comparison_file_path, 'w') as f:
            f.write(comparison_latex)
        print(f"  - Success: Saved comparison report to '{comparison_file_path}'")
        
    except Exception as e:
        print(f"  [Error] Failed to write comparison file: {e}", file=sys.stderr)

    # --- MODIFIED: This loop now saves files using the file_prefix to match your \input ---
    # --- 2. Generate and save a separate .tex file for each algorithm ---
    # We iterate over the MAP to get the original file_prefix for the filename
    for file_prefix, algo_name in ALGO_NAME_MAP.items():
        if algo_name not in all_algo_data:
            print(f"  [Info] Skipping detail file for '{algo_name}' (no data was parsed).")
            continue
        
        data = all_algo_data[algo_name]
        try:
            # Generate the table fragment using the "pretty name"
            detail_latex = generate_algo_table_fragment(algo_name, data)
            
            # Create a safe filename FROM THE FILE_PREFIX to match \input
            safe_filename = "detail_" + re.sub(r'[^a-zA-Z0-9]', '', file_prefix) + ".tex"
            detail_file_path = OUTPUT_FOLDER / safe_filename
            
            with open(detail_file_path, 'w') as f:
                f.write(detail_latex)
            print(f"  - Success: Saved detail report to '{detail_file_path}'")
            
        except Exception as e:
            print(f"  [Error] Failed to write detail file for '{algo_name}': {e}", file=sys.stderr)
    # --- END OF MODIFIED LOOP ---

    print("\n--- Generation Complete! ---")
    print("You can now upload the .tex files from")
    print(f"'{OUTPUT_FOLDER}'")
    print("to your Overleaf project and use \input{...} to include them.")
    print("Remember to add \\usepackage{longtable} to your main .tex file!")

if __name__ == "__main__":
    main()