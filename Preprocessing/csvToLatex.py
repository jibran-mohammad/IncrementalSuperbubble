import pandas as pd
import sys
import os

def generate_latex_longtable(csv_filename, output_filename="haplotype_longtable.tex"):
    """
    Reads a CSV, selects specific columns, and generates a multi-page compatible 
    LaTeX 'longtable' file, ensuring vertical column lines are included.

    Args:
        csv_filename (str): The name of the input CSV file ("haplotype_build_comparison.csv").
        output_filename (str): The name of the output LaTeX (.tex) file.
    """
    try:
        # --- 1. Load the data ---
        df = pd.read_csv(csv_filename)

        # --- 2. Select and validate columns ---
        required_columns = [
            'Step', 
            'HaplotypeName', 
            'Nodes', 
            'Edges', 
            'SuperbubbleCount'
        ]
        
        missing_cols = [col for col in required_columns if col not in df.columns]
        if missing_cols:
            print(f"Error: Missing required columns: {', '.join(missing_cols)}")
            sys.exit(1)

        df_selected = df[required_columns]
        num_cols = len(required_columns)
        
        # Column format: |c|c|c|c|c| ensures vertical lines between all columns
        column_format = '|' + '|'.join(['c'] * num_cols) + '|'

        # --- 3. Generate the LaTeX content using 'longtable' ---
        
        latex_content = []

        # Start of the longtable environment
        # Note: 'longtable' cannot be placed inside a 'table' float environment
        latex_content.append(r"\begin{longtable}{" + column_format + "}")
        
        # Caption and Header for the FIRST page
        latex_content.append(r"    \caption{Haplotype Build Comparison Data} \\")
        latex_content.append(r"    \label{tab:haplotype_data}")
        latex_content.append(r"    \hline")
        
        # Header row - Escape underscores and use \textbf{}
        header_row = [col.replace('_', r'\_') for col in df_selected.columns]
        latex_content.append("    " + " & ".join([f"\\textbf{{{h}}}" for h in header_row]) + r" \\ \hline")
        
        # Defines the end of the first page header
        latex_content.append(r"    \endfirsthead")
        
        # Caption and Header for SUBSEQUENT pages (repeated on page break)
        latex_content.append(r"    \multicolumn{" + str(num_cols) + r"}{|c|}{\textbf{\tablename\ \thetable\ -- Continued}} \\")
        latex_content.append(r"    \hline")
        latex_content.append("    " + " & ".join([f"\\textbf{{{h}}}" for h in header_row]) + r" \\ \hline")
        latex_content.append(r"    \endhead")
        
        # Footer when the table continues
        latex_content.append(r"    \hline")
        latex_content.append(r"    \multicolumn{" + str(num_cols) + r"}{|r|}{{Continued on next page}} \\")
        latex_content.append(r"    \endfoot")
        
        # Footer for the last page
        latex_content.append(r"    \hline")
        latex_content.append(r"    \endlastfoot")

        # Data rows
        for index, row in df_selected.iterrows():
            # Format the data, ensuring HaplotypeName's underscores are escaped for LaTeX
            row_data = [str(cell).replace('_', r'\_') if col == 'HaplotypeName' else str(cell) 
                        for col, cell in zip(df_selected.columns, row)]
            latex_content.append("    " + " & ".join(row_data) + r" \\ \hline")

        # End of the longtable environment
        latex_content.append(r"\end{longtable}")

        # --- 4. Write to the output file ---
        with open(output_filename, 'w') as f:
            f.write('\n'.join(latex_content))

        print(f"✅ Successfully generated multi-page LaTeX table in file: **{output_filename}**")
        print("Note: Remember to use \\usepackage{longtable} in your main LaTeX document.")

    except FileNotFoundError:
        print(f"Error: Input file '{csv_filename}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)

# --- Execution Block ---
if __name__ == "__main__":
    csv_file_name = "/home/dell/Research/Superbubble/Experiment/IncSuperbubble/output/Ecoli/timeFile/haplotype_build_comparison.csv"
    
    if os.path.exists(csv_file_name):
        generate_latex_longtable(csv_file_name)
    else:
        print(f"Error: The file **{csv_file_name}** was not found in the current directory.")