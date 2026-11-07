"""
split_fasta.py

Splits a multi-chromosome reference FASTA (e.g. T2T-CHM13, CHM13Y.fa) into
one file per chromosome, in the layout the graph pipeline expects
(scripts/graphRunner.py -d Human):

    <outdir>/chromosome_<c>/chromosome_<c>.fa

The chromosome name comes from the FASTA header: ">chr1 ..." -> 1,
">chrX" -> X, ">chrMT" / ">MT" -> M.

=====================================================================
HOW TO RUN
=====================================================================

    # writes input/Human/chromosome_1/chromosome_1.fa, ... (next to the FASTA)
    python3 Preprocessing/split_fasta.py input/Human/CHM13Y.fa

    # another output folder; gzipped input also works
    python3 Preprocessing/split_fasta.py CHM13Y.fa.gz --outdir input/Human

Options:
    fasta       the reference FASTA (.fa/.fna/.fasta, optionally .gz)
    --outdir    folder for the chromosome_<c>/ folders
                (default: the folder of the FASTA)
"""
import argparse
import gzip
import os
import sys


def main():
    ap = argparse.ArgumentParser(description="Split a reference FASTA into one file per chromosome")
    ap.add_argument("fasta", help="reference FASTA (optionally .gz)")
    ap.add_argument("--outdir", help="folder for the chromosome_<c>/ folders (default: the folder of the FASTA)")
    args = ap.parse_args()

    input_fasta = args.fasta
    if not os.path.isfile(input_fasta):
        sys.exit(f"❌ Error: the file {input_fasta} was not found.")
    base_dir = args.outdir or os.path.dirname(os.path.abspath(input_fasta))
    opener = gzip.open if input_fasta.endswith(".gz") else open

    print(f"Starting to split {input_fasta}...")
    output_file_handle = None
    with opener(input_fasta, "rt") as f_in:
        for line in f_in:
            if line.startswith(">"):
                if output_file_handle:
                    output_file_handle.close()

                raw_id = line.strip()[1:].split()[0]
                base_id = raw_id[3:] if raw_id.startswith("chr") else raw_id
                if base_id == "MT":
                    base_id = "M"

                # chromosome_<c>/chromosome_<c>.fa
                dir_name_base = "chromosome_" + base_id
                full_dir_path = os.path.join(base_dir, dir_name_base)
                os.makedirs(full_dir_path, exist_ok=True)
                output_path = os.path.join(full_dir_path, dir_name_base + ".fa")
                output_file_handle = open(output_path, "w")
                print(f"  -> Writing to {output_path}")

            if output_file_handle:
                output_file_handle.write(line)

    if output_file_handle:
        output_file_handle.close()
    print(f"\n✅ Success! Finished splitting {input_fasta}.")


if __name__ == "__main__":
    main()
