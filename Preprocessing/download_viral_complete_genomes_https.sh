#!/usr/bin/env bash
# Usage: ./download_viral_complete_genomes_https.sh assemblies.txt out_dir
# assemblies.txt: one NCBI assembly *folder* URL per line (ftp:// or https://, with or without trailing slash)
#
# For each assembly:
#   1. Check *_assembly_report.txt for "Assembly level: Complete Genome"
#   2. If complete, download *_genomic.fna(.gz) files.

set -euo pipefail

LIST=${1:-viral_list.txt}
OUT=${2:-viral_genomes_out}
mkdir -p "$OUT"

to_https() {
  # ftp://…  -> https://…
  sed -E 's#^ftp://#https://#'
}

while IFS= read -r base || [ -n "$base" ]; do
  [ -z "$base" ] && continue
  base="${base%/}/"
  base_https="$(printf "%s" "$base" | to_https)"
  asm="$(basename "$base_https")"
  dest="$OUT/$asm"

  # If already downloaded, skip
  if [ -d "$dest" ]; then
    echo "SKIPPING $asm: Directory '$dest' already exists."
    continue
  fi

  # --- 1. Check assembly_report for 'Complete Genome' ---
  report_url="${base_https}${asm}_assembly_report.txt"
  echo "Checking completeness for $asm ..."

  if ! curl -fsSL "$report_url" | grep -qi 'Assembly level:[[:space:]]*Complete Genome'; then
    echo "❌ Skipping $asm: not a 'Complete Genome' assembly (or no report)."
    continue
  fi

  # --- 2. Find *_genomic.fna(.gz) in the base directory ---
  echo "  → $asm is 'Complete Genome'; listing genomic FASTA files..."
  files=$(curl -fsSL "$base_https" \
          | grep -Eo '[A-Za-z0-9._-]+_genomic\.fna(\.gz)?' \
          | sort -u)

  if [ -z "$files" ]; then
    echo "⚠️  $asm is complete but no *_genomic.fna(.gz) files found at $base_https"
    continue
  fi

  mkdir -p "$dest"
  echo "⬇️  $asm : $(echo "$files" | wc -l) genomic FASTA file(s)"

  while read -r f; do
    [ -z "$f" ] && continue
    echo "  → Downloading $f"
    curl -fsSLO --output-dir "$dest" "${base_https}${f}"
  done <<< "$files"

done < "$LIST"

echo "✅ Done. Output in: $OUT"
