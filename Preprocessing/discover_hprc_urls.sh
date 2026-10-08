#!/usr/bin/env bash
# discover_hprc_urls.sh  - produce a list of NCBI FTP assembly folder URLs for HPRC entries
# Usage: ./discover_hprc_urls.sh -n 3 -o list_test.txt
set -euo pipefail

N=3
OUT="hprc_list.txt"
HUB_ROOT="https://hgdownload.soe.ucsc.edu/hubs/HPRC/"

while getopts "n:o:h" opt; do
  case $opt in
    n) N="$OPTARG" ;;
    o) OUT="$OPTARG" ;;
    h) echo "Usage: $0 [-n NUM] [-o out_file]"; exit 0 ;;
    *) echo "Unknown option"; exit 1 ;;
  esac
done

tmp=$(mktemp)
trap 'rm -f "$tmp"' EXIT

echo "Fetching hub index from $HUB_ROOT ..."
curl -sSL "$HUB_ROOT" -o "$tmp"

# extract tokens like: GCA_044165215.1_HG00097_hap1_hprc_f2
grep -Eo 'GCA_[0-9]+\.[0-9]+_[A-Za-z0-9._-]+' "$tmp" \
  | awk '!seen[$0]++' > "$tmp.tokens"

if [ ! -s "$tmp.tokens" ]; then
  echo "No GCA tokens found on hub page; inspect $tmp"
  exit 1
fi

# produce FTP paths
: > "$OUT.tmp"
count=0
while read -r token; do
  # optionally limit
  count=$((count+1))
  if [ "$N" -gt 0 ] && [ "$count" -gt "$N" ]; then break; fi

  # accession part (GCA_XXXXXXXXX.Y)
  acc=$(printf "%s" "$token" | sed -E 's/^((GCA_[0-9]+\.[0-9]+)).*/\1/')
  # numeric portion (only digits)
  num=$(printf "%s" "$acc" | sed -E 's/GCA_([0-9]+)\..*/\1/')

  # pad left with zeros until length is multiple of 3
  while [ $(( ${#num} % 3 )) -ne 0 ]; do num="0${num}"; done

  d1=${num:0:3}; d2=${num:3:3}; d3=${num:6:3}
  printf "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/%s/%s/%s/%s/\n" "$d1" "$d2" "$d3" "$token" >> "$OUT.tmp"
done < "$tmp.tokens"

mv "$OUT.tmp" "$OUT"
echo "Wrote $(wc -l < "$OUT") URLs to $OUT"
sed -n '1,20p' "$OUT"
