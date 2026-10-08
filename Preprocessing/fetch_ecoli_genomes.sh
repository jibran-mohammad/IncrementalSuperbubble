#!/usr/bin/env bash
set -euo pipefail

if [ -z "$1" ]; then
    echo "Error: Please provide a base directory."
    echo "Usage: $0 /path/to/your/project_folder"
    exit 1
fi

# Use realpath to get the full, absolute path
BASE_DIR=$(realpath "$1")
echo "Target base directory set to: $BASE_DIR"

# --- Original script logic, now using BASE_DIR ---
WORKDIR="$BASE_DIR/ecoli_100"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

ASM_SUM="assembly_summary_refseq.txt"
OUT_CAND="ecoli_candidate_ftps.txt"
OUT_URLS="100_urls.txt"
EXCLUDE_REF="GCF_000005845.2"

echo "1) Downloading RefSeq assembly summary..."
wget -q -O "$ASM_SUM" "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"

echo "2) Filtering for Escherichia coli + Complete Genome..."
header_line=$(grep -m1 -n "^# assembly_accession" -n "$ASM_SUM" | cut -d: -f1 || true)
if [ -z "$header_line" ]; then
  header_line=1
fi

awk -F'\t' '
  BEGIN{OFS="\t"}
  /^#/ { next }
  {
    acc=$1;
    asm_level=$12;
    orgname=$8;
    ftp=$20;
    if(ftp=="") ftp=$(NF);
    print acc,orgname,asm_level,ftp
  }' "$ASM_SUM" > asm_summary_parsed.tsv

awk -F'\t' 'tolower($2) ~ /escherichia coli/ && tolower($3)=="complete genome" {print $0}' asm_summary_parsed.tsv \
  | awk -v excl="$EXCLUDE_REF" -F'\t' '$1 != excl {print $0}' > "$OUT_CAND"

wc -l "$OUT_CAND" || true
echo "Found candidate lines (complete genome Escherichia coli)."

echo "3) Build genomic.fna.gz URLs and verify existence..."
> "$OUT_URLS"

python3 - <<'PY'
import sys,os,urllib.request,urllib.error
cand_file="ecoli_candidate_ftps.txt"
out_file="100_urls.txt"
seen_strains=set()
count=0
max_needed=3
def exists(url):
    try:
        req=urllib.request.Request(url, method='HEAD')
        with urllib.request.urlopen(req, timeout=15) as r:
            return r.status==200
    except Exception:
        return False

with open(cand_file) as fh, open(out_file,'w') as out:
    for line in fh:
        acc,org,asm_level,ftp=line.rstrip("\n").split("\t")
        asm_name = ftp.rstrip("/").split("/")[-1]
        fna_name = f"{asm_name}_genomic.fna.gz"
        url = ftp.rstrip("/") + "/" + fna_name
        strain_key=org.lower().strip()
        if strain_key in seen_strains:
            continue
        if exists(url):
            out.write(f"{acc}\t{org}\t{url}\n")
            seen_strains.add(strain_key)
            count += 1
            print(f"[{count}] OK: {acc} {org}")
            sys.stdout.flush()
            if count >= max_needed:
                break
        else:
            print(f"MISS: {acc} -- {url}", file=sys.stderr)
PY

echo "Done. '$OUT_URLS' contains up to 100 verified URLs."
echo "Preview (first 10):"
head -n 10 "$OUT_URLS" || true


echo "4) Creating genomes download directory..."
# Create the 'genomes' directory *inside* the BASE_DIR
GENOME_DIR="$BASE_DIR/genomes"
mkdir -p "$GENOME_DIR"

echo "5) Downloading all .fna.gz files into $GENOME_DIR..."
# Change to the genomes directory before downloading
cd "$GENOME_DIR"

# Get the full path to the URL list
URL_FILE_PATH="$WORKDIR/$OUT_URLS"

# Use your command to pipe the 3rd column (URLs) to wget
cut -f3 "$URL_FILE_PATH" | wget -c -i -

echo "6) Unzipping all downloaded .fna.gz files..."
# We are still in the $GENOME_DIR
gunzip *.fna.gz

echo "---"
echo "✅ All steps complete."
echo "URL list: $URL_FILE_PATH"
echo "Genomes:  $GENOME_DIR"