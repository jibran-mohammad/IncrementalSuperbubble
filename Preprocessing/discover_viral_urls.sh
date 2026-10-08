#!/usr/bin/env bash
# discover_viral_urls.sh  - produce a list of NCBI FTP assembly folder URLs for viral entries
# Usage:
#   ./discover_viral_urls.sh -n 10 -o viral_list.txt
#   ./discover_viral_urls.sh -n 0  -o viral_list_all.txt   # 0 = no limit

set -euo pipefail

# Number of assemblies to output (0 = unlimited)
N=3
OUT="viral_list.txt"
HUB_ROOT="https://hgdownload.soe.ucsc.edu/hubs/viral/"

while getopts "n:o:h" opt; do
  case $opt in
    n) N="$OPTARG" ;;
    o) OUT="$OPTARG" ;;
    h) echo "Usage: $0 [-n NUM] [-o out_file]"; exit 0 ;;
    *) echo "Unknown option"; exit 1 ;;
  esac
done

tmp=$(mktemp)
trap 'rm -f "$tmp" "$tmp.tokens" 2>/dev/null || true' EXIT

echo "Fetching viral hub index from $HUB_ROOT ..."
curl -sSL "$HUB_ROOT" -o "$tmp"

# Extract tokens like:
#   GCA_044165215.1_SomeVirus_something
#   GCF_012345678.1_SomeVirus_something
grep -Eo 'GC[AF]_[0-9]+\.[0-9]+_[A-Za-z0-9._-]+' "$tmp" \
  | awk '!seen[$0]++' > "$tmp.tokens"

if [ ! -s "$tmp.tokens" ]; then
  echo "No GC[AF]_ tokens found on hub page; inspect $tmp"
  exit 1
fi

: > "$OUT.tmp"
count=0
while read -r token; do
  [ -z "$token" ] && continue
  count=$((count+1))
  # N == 0  → no upper limit
  if [ "$N" -gt 0 ] && [ "$count" -gt "$N" ]; then
    break
  fi

  # accession part (GCA_XXXXXXXXX.Y or GCF_XXXXXXXXX.Y)
  acc=$(printf "%s" "$token" | sed -E 's/^((GC[AF]_[0-9]+\.[0-9]+)).*/\1/')
  # numeric portion (only digits)
  num=$(printf "%s" "$acc" | sed -E 's/GC[AF]_([0-9]+)\..*/\1/')

  # pad left with zeros until length is multiple of 3
  while [ $(( ${#num} % 3 )) -ne 0 ]; do num="0${num}"; done

  d1=${num:0:3}; d2=${num:3:3}; d3=${num:6:3}

  printf "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/%s/%s/%s/%s/%s/\n" \
    "${acc%_*}" "$d1" "$d2" "$d3" "$token" >> "$OUT.tmp"

  # NOTE:
  # The path pattern above assumes the NCBI directory layout used for GCA/GCF.
  # If you see 404s, run with a small -n and adjust the printf pattern
  # based on an example URL from NCBI.
done < "$tmp.tokens"

mv "$OUT.tmp" "$OUT"
echo "Wrote $(wc -l < "$OUT") URLs to $OUT"
sed -n '1,20p' "$OUT"
