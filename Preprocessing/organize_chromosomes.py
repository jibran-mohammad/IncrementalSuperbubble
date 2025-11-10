#!/usr/bin/env python3

"""
Organize Haplotype Assemblies by Chromosome

This script reorganizes a directory of haplotype assemblies into a 
new directory structure sorted by chromosome.

Expected Input Structure:
--------------------------
source_directory/
    ├── haplotype_1_dir/
    │   ├── chr1.fasta
    │   ├── chr2.fasta
    │   └── chrM.fasta
    ├── haplotype_2_dir/
    │   ├── chr1.fasta
    │   ├── chrX.fasta
    │   └── chrY.fasta
    └── ...

Resulting Output Structure:
---------------------------
output_directory/
    ├── chromosome_1/
    │   ├── haplotype_1_dir.fasta
    │   └── haplotype_2_dir.fasta
    ├── chromosome_2/
    │   └── haplotype_1_dir.fasta
    ├── chromosome_X/
    │   └── haplotype_2_dir.fasta
    ├── chromosome_Y/
    │   └── haplotype_2_dir.fasta
    ├── chromosome_M/
    │   └── haplotype_1_dir.fasta
    └── ...

The script identifies chromosomes (1-22, X, Y, M) from the fasta filenames
and uses the haplotype directory name for the new, copied fasta file.

Usage:
------

1. Run from your terminal:

   # Basic usage (output will be in 'organized_chromosomes')
   python3 organize_chromosomes.py /path/to/your/input_directory

   # Specify a custom output directory
   python3 organize_chromosomes.py /path/to/your/input_directory -o /path/to/your/output_directory

   # Get help
   python3 organize_chromosomes.py -h
"""

import os
import sys
import shutil
import argparse
import re
import gzip  # Import the gzip module
from typing import Optional

def get_standard_chromosome_dir(filename: str) -> Optional[str]:
    """
    Normalizes a fasta filename to a standard chromosome directory name.
    
    Examples:
        "chr1.fasta" -> "chromosome_1"
        "01.fa"      -> "chromosome_1"
        "Y.fa"       -> "chromosome_Y"
        "chrM.fna"   -> "chromosome_M"
        "mito.fa"    -> "chromosome_M"
        "assembly.chrMT.fasta" -> "chromosome_M"

    Args:
        filename: The original fasta filename (e.g., "chr1.fasta").

    Returns:
        A standardized directory name (e.g., "chromosome_1") or None
        if no chromosome pattern is recognized.
    """
    # Get the part before the first dot (e.g., "chr1.fasta.gz" -> "chr1")
    basename = filename.split('.')[0]
    
    # Regex to find the chromosome identifier at the end of the basename
    # Matches: 
    # - Optional 'chr' prefix
    # - 1-22 (with optional leading 0, e.g., 01, 1, 22)
    # - X, Y, M, or MT
    # - Anchored to the end of the basename ($)
    match = re.search(r'(?:chr)?(2[0-2]|1[0-9]|0?[1-9]|X|Y|M|MT)$', basename, re.IGNORECASE)
    
    chrom_identifier = None
    if match:
        # Get the captured part (e.g., "01", "1", "X", "MT")
        chrom_identifier = match.group(1).upper() 
    elif 'mitochondrion' in basename.lower() or 'mito' in basename.lower():
        # Handle cases where the name is "mitochondrion.fasta"
        chrom_identifier = 'M'
    
    if chrom_identifier:
        # Standardize 'MT' to 'M'
        if chrom_identifier == 'MT':
            chrom_identifier = 'M'
        
        # Standardize '01', '02' etc. to '1', '2'
        # while keeping '10', '20' etc. as is.
        if chrom_identifier.isdigit():
            chrom_identifier = str(int(chrom_identifier))
            
        return f"chromosome_{chrom_identifier}"
    
    # If no pattern matches, return None
    return None

def main():
    """
    Main function to parse arguments and run the organization process.
    """
    parser = argparse.ArgumentParser(
        description="Organize haplotype assemblies by chromosome.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__ # Use the script's docstring as the help epilog
    )
    
    parser.add_argument(
        "input_dir",
        type=str,
        help="The source directory containing haplotype subdirectories (e.g., 'chromosomes_test_out')."
    )
    
    parser.add_argument(
        "-o", "--output_dir",
        type=str,
        default="organized_chromosomes",
        help="The destination directory to create the organized structure. (Default: 'organized_chromosomes')"
    )
    
    args = parser.parse_args()
    
    input_dir = args.input_dir
    output_dir = args.output_dir

    # --- 1. Validate Input Directory ---
    if not os.path.isdir(input_dir):
        print(f"Error: Input directory not found or is not a directory: {input_dir}", file=sys.stderr)
        sys.exit(1)
        
    # --- 2. Create Output Directory ---
    try:
        os.makedirs(output_dir, exist_ok=True)
        print(f"Organizing chromosomes into: {os.path.abspath(output_dir)}")
    except OSError as e:
        print(f"Error: Could not create output directory: {output_dir}\n{e}", file=sys.stderr)
        sys.exit(1)

    # --- 3. Find Haplotype Directories ---
    try:
        haplotype_dirs = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]
    except OSError as e:
        print(f"Error: Could not read input directory {input_dir}: {e}", file=sys.stderr)
        sys.exit(1)

    if not haplotype_dirs:
        print(f"Warning: No subdirectories (haplotypes) found in {input_dir}.")
        return

    print(f"Found {len(haplotype_dirs)} potential haplotype directories.")
    
    # --- 4. Process Each Haplotype Directory ---
    file_count = 0
    skipped_files = 0
    
    for hap_dir_name in haplotype_dirs:
        hap_dir_path = os.path.join(input_dir, hap_dir_name)
        print(f"\nProcessing haplotype: {hap_dir_name}...")
        
        try:
            # Find all common fasta extensions, including gzipped files
            all_files = os.listdir(hap_dir_path)
            fasta_files = []
            for f in all_files:
                f_lower = f.lower()
                is_gz = f_lower.endswith('.gz')
                is_fasta = f_lower.endswith(('.fa', '.fasta', '.fna'))
                
                if is_gz:
                    # Check for .fa.gz, .fasta.gz, .fna.gz
                    if (f_lower.endswith('.fa.gz') or 
                        f_lower.endswith('.fasta.gz') or 
                        f_lower.endswith('.fna.gz')):
                        fasta_files.append(f)
                elif is_fasta:
                    # Handle non-gzipped fasta files
                    fasta_files.append(f)

        except OSError as e:
            print(f"  Warning: Could not read directory {hap_dir_path}. Skipping. Error: {e}", file=sys.stderr)
            continue

        if not fasta_files:
            print(f"  Info: No fasta files (.fa, .fasta, .fna) found in {hap_dir_path}.")
            continue

        # --- 5. Process Each Fasta File ---
        for filename in fasta_files:
            source_filepath = os.path.join(hap_dir_path, filename)
            
            # Get the standard destination directory name (e.g., "chromosome_1")
            std_chrom_dir_name = get_standard_chromosome_dir(filename)
            
            if not std_chrom_dir_name:
                print(f"  Warning: Skipping file. Could not normalize chromosome name for: {filename}")
                skipped_files += 1
                continue
            
            # --- 6. Create Destination Directory and File Paths ---
            
            # Full path for the chromosome-specific directory (e.g., ./organized_chromosomes/chromosome_1)
            dest_chrom_dir_path = os.path.join(output_dir, std_chrom_dir_name)
            os.makedirs(dest_chrom_dir_path, exist_ok=True)
            
            # New filename will be the haplotype directory name + .fasta
            # (e.g., GCA_044164745.1_HG00097_hap2_hprc_f2.fasta)
            dest_filename = f"{hap_dir_name}.fasta"
            dest_filepath = os.path.join(dest_chrom_dir_path, dest_filename)
            
            if os.path.exists(dest_filepath):
                print(f"  Info: Overwriting existing file: {os.path.join(std_chrom_dir_name, dest_filename)}")
            
            # --- 7. Copy and Unzip the File (if needed) ---
            try:
                # Check if the source file is gzipped
                if source_filepath.lower().endswith('.gz'):
                    # Open the gzipped file for reading in binary mode
                    with gzip.open(source_filepath, 'rb') as f_in:
                        # Open the destination file for writing in binary mode
                        with open(dest_filepath, 'wb') as f_out:
                            # Copy and decompress data
                            shutil.copyfileobj(f_in, f_out)
                    print(f"  [OK] {filename:20} -> {os.path.join(std_chrom_dir_name, dest_filename)} (unzipped)")
                else:
                    # If not gzipped, just copy it as before
                    shutil.copy2(source_filepath, dest_filepath)
                    print(f"  [OK] {filename:20} -> {os.path.join(std_chrom_dir_name, dest_filename)}")
                
                file_count += 1
            except shutil.SameFileError:
                 print(f"  Warning: Source and destination are the same file: {source_filepath}", file=sys.stderr)
            except (OSError, gzip.BadGzipFile) as e:
                print(f"  Error: Could not copy or unzip {filename}. Error: {e}", file=sys.stderr)

    print(f"\n---------------------------------")
    print(f"Organization complete.")
    print(f"  Total files copied: {file_count}")
    print(f"  Total files skipped: {skipped_files}")
    print(f"---------------------------------")

if __name__ == "__main__":
    main()


