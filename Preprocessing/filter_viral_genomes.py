#!/usr/bin/env python3
"""
Filter viral assembly directories:

- Remove files whose names contain: cds, rna, protein
- Extract *_genomic.fna.gz → *_genomic.fna and delete .gz
- Collect only single-sequence *_genomic.fna into: cleanGenomes/
  (multi-segment/multi-contig genomes are SKIPPED)

Usage:
    python filter_extract_collect_viral_genomes.py /path/to/viral_genomes
"""

import os
import sys
import gzip
import shutil

PATTERNS_TO_DELETE = ("cds", "rna", "protein")


def delete_unwanted_files(root_dir: str) -> None:
    """Delete files containing cds/rna/protein in their names (case-insensitive)."""
    for dirpath, _, filenames in os.walk(root_dir):
        for fname in filenames:
            lower = fname.lower()
            if any(p in lower for p in PATTERNS_TO_DELETE):
                path = os.path.join(dirpath, fname)
                print(f"Removing: {path}")
                try:
                    os.remove(path)
                except Exception as e:
                    print(f"  [WARN] Could not remove {path}: {e}")


def extract_genomic_gz(root_dir: str) -> None:
    """Extract *_genomic.fna.gz and delete original gzip."""
    for dirpath, _, filenames in os.walk(root_dir):
        for fname in filenames:
            if fname.lower().endswith("genomic.fna.gz"):
                gz_path = os.path.join(dirpath, fname)
                out_path = gz_path[:-3]

                print(f"Extracting: {gz_path} -> {out_path}")
                try:
                    with gzip.open(gz_path, "rb") as f_in, open(out_path, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)
                    os.remove(gz_path)
                except Exception as e:
                    print(f"  [WARN] Failed unzip {gz_path}: {e}")


def is_single_sequence_fasta(path: str) -> bool:
    """Return True if FASTA file has exactly one '>' header; False otherwise."""
    headers = 0
    try:
        with open(path, "r") as f:
            for line in f:
                if line.startswith(">"):
                    headers += 1
                    if headers > 1:
                        return False
    except Exception as e:
        print(f"  [WARN] Could not read {path}: {e}")
        return False
    return headers == 1


def collect_genomes(root_dir: str) -> None:
    """
    Collect all single-sequence *_genomic.fna into cleanGenomes/.
    Multi-segment (multi-header) files are skipped.
    """
    output_dir = os.path.join(root_dir, "cleanGenomes")
    os.makedirs(output_dir, exist_ok=True)

    kept = 0
    skipped = 0

    for dirpath, _, filenames in os.walk(root_dir):
        for fname in filenames:
            if fname.lower().endswith("genomic.fna"):
                full_path = os.path.join(dirpath, fname)

                if not is_single_sequence_fasta(full_path):
                    print(f"Skipping multi-segment genome: {full_path}")
                    skipped += 1
                    # If you also want to DELETE them, uncomment:
                    # os.remove(full_path)
                    continue

                dirname = os.path.basename(dirpath)
                out_name = f"{dirname}.fna"
                out_path = os.path.join(output_dir, out_name)

                print(f"Collecting: {full_path} → {out_path}")
                try:
                    shutil.copy(full_path, out_path)
                    kept += 1
                except Exception as e:
                    print(f"  [WARN] Could not copy {full_path}: {e}")

    print(f"\n📦 Collected {kept} single-contig genomes into: {output_dir}")
    print(f"🚫 Skipped {skipped} multi-segment genomes.\n")


def main():
    if len(sys.argv) != 2:
        print("Usage: python filter_extract_collect_viral_genomes.py /path/to/viral_genomes")
        sys.exit(1)

    root_dir = sys.argv[1]
    if not os.path.isdir(root_dir):
        print("ERROR: Provided directory does not exist")
        sys.exit(1)

    print(f"Processing viral genomes under: {root_dir}\n")

    delete_unwanted_files(root_dir)
    extract_genomic_gz(root_dir)
    collect_genomes(root_dir)

    print("✨ Done. cleanGenomes/ now contains only single-sequence genomes.\n")


if __name__ == "__main__":
    main()
