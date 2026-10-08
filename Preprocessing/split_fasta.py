import os
import sys

# --- Configuration ---
# Make sure this path is correct, even if it's in another folder
# e.g., "data/CHM13Y.fa" or "/home/user/genomes/CHM13Y.fa"
input_fasta = "/home/dell/Research/Superbubble/Experiment/IncSuperbubble/input/Human1/CHM13Y.fa"
# ---------------------

# <-- NEW LINE: Get the directory where the FASTA file is located.
base_dir = os.path.dirname(input_fasta) 

output_file_handle = None

print(f"Starting to split {input_fasta}...")

try:
    with open(input_fasta, 'r') as f_in:
        for line in f_in:
            if line.startswith('>'):
                if output_file_handle:
                    output_file_handle.close()
                
                raw_id = line.strip()[1:].split()[0] 
                
                base_id = ""
                if raw_id.startswith('chr'):
                    base_id = raw_id[3:]
                else:
                    base_id = raw_id
                
                if base_id == 'MT':
                    base_id = 'M'
                
                # Create the base name for dir and file
                dir_name_base = "chromosome_" + base_id
                file_name = dir_name_base + ".fa"
                
                # <-- UPDATED LINE: Join the base_dir with the new chromosome dir
                full_dir_path = os.path.join(base_dir, dir_name_base)
                
                # 5. Create the directory using the full path
                os.makedirs(full_dir_path, exist_ok=True) # <-- UPDATED LINE
                
                # 6. Define the full output path
                output_path = os.path.join(full_dir_path, file_name)
                output_file_handle = open(output_path, 'w')
                
                print(f"  -> Writing to {output_path}")

            if output_file_handle:
                output_file_handle.write(line)

    if output_file_handle:
        output_file_handle.close()

    print(f"\n✅ Success! Finished splitting {input_fasta}.")

except FileNotFoundError:
    print(f"❌ Error: The file {input_fasta} was not found.")
    print("Please make sure the path in the 'input_fasta' variable is correct.")
except Exception as e:
    print(f"An unexpected error occurred: {e}")