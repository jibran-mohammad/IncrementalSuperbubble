#!/usr/bin/env bash
# Usage: ./download_chromosomes_https.sh assemblies.txt out_dir
# assemblies.txt: one NCBI assembly *folder* URL per line (ftp:// or https://, with or without trailing slash)

set -euo pipefail
LIST=${1:-list_for_test.txt}
OUT=${2:-chromosomes_out}
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

  chr_url="${base_https}${asm}_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/"

  # Check that folder exists (HTTP 200)
  if ! curl -fsSI "$chr_url" >/dev/null; then
    echo "❌ No chromosome FASTA folder for $asm"
    continue
  fi

  # Parse filenames (robust; handles plain listings)
  files=$(curl -fsSL "$chr_url" \
          | grep -Eo 'chr(1?[0-9]|2[0-2]|X|Y|MT)\.fna(\.gz)?' \
          | sort -u)

  if [ -z "$files" ]; then
    echo "⚠️  Found folder but no chr*.fna(.gz) for $asm"
    continue
  fi

  dest="$OUT/$asm"
  mkdir -p "$dest"
  echo "⬇️  $asm : $(echo "$files" | wc -l) chromosomes"
  while read -r f; do
    [ -z "$f" ] && continue
    echo "  → $f"
    curl -fsSLO --output-dir "$dest" "$chr_url$f"
  done <<< "$files"
done < "$LIST"

echo "✅ Done. Output in: $OUT"
