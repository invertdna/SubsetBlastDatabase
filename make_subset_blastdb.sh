#!/usr/bin/env bash
# make_subset_blastdb.sh
# Fetch accession numbers from NCBI nuccore and build a standalone subset BLAST database.
# Bash equivalent of make_subset_blastdb.R — requires only BLAST+ and edirect (no R).
#
# Usage:
#   bash make_subset_blastdb.sh <blastdb> <taxon> <query> <title> [output_dir]
#
# Arguments:
#   blastdb    : path (with db prefix) to an existing local BLAST db
#                  e.g. /Volumes/Clupea/core_nt/core_nt
#   taxon      : short name used for output file naming — no spaces (e.g. Sebastes)
#   query      : NCBI esearch query string
#   title      : human-readable title embedded in the output database metadata
#   output_dir : (optional) output directory; defaults to the taxon value

set -euo pipefail

CHUNK_SIZE=100000

# ---- parse args --------------------------------------------------------
if [[ $# -lt 4 ]]; then
  cat >&2 <<EOF
Usage: bash make_subset_blastdb.sh <blastdb> <taxon> <query> <title> [output_dir]
  blastdb    path to existing local BLAST db
  taxon      short name for file naming (no spaces)
  query      NCBI esearch query string
  title      title for the output database
  output_dir (optional) output directory; default = taxon
EOF
  exit 1
fi

db_path="$1"
taxon="$2"
query="$3"
title="$4"
out_dir="${5:-$taxon}"

# ---- locate subset_blastdb_by_accession.sh next to this script ---------
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
subset_script="${script_dir}/subset_blastdb_by_accession.sh"
if [[ ! -f "$subset_script" ]]; then
  echo "ERROR: Cannot find subset_blastdb_by_accession.sh in: ${script_dir}" >&2
  exit 1
fi

# ---- ensure edirect is on PATH -----------------------------------------
export PATH="${PATH}:${HOME}/edirect"

# =======================================================================
# STEP 1 — esearch
# =======================================================================
echo ""
echo "--- Step 1: Querying NCBI nuccore ---"
echo "  query: ${query}"

search_xml="$(mktemp /tmp/esearch_XXXXXX.xml)"
trap 'rm -f "$search_xml"' EXIT

echo "  cmd: esearch -db nuccore -query '${query}'"
esearch -db nuccore -query "$query" > "$search_xml"

total=$(grep -o '<Count>[0-9]*</Count>' "$search_xml" | head -1 | grep -o '[0-9]*' || true)
if [[ -z "$total" || "$total" -eq 0 ]]; then
  echo "ERROR: No records found or could not parse count from esearch output." >&2
  cat "$search_xml" >&2
  exit 1
fi
echo "  Total matching records: ${total}"

# =======================================================================
# STEP 2 — efetch accessions in chunks
# =======================================================================
echo ""
echo "--- Step 2: Fetching accessions ---"

acc_dir="${out_dir}/accessions"
mkdir -p "$acc_dir"

n_chunks=$(( (total + CHUNK_SIZE - 1) / CHUNK_SIZE ))
chunk_num=1
start=0
chunk_files=()

while [[ "$start" -lt "$total" ]]; do
  end=$(( start + CHUNK_SIZE - 1 ))
  if [[ "$end" -ge "$total" ]]; then
    end=$(( total - 1 ))
  fi

  echo "  Chunk ${chunk_num}/${n_chunks}  (records ${start} – ${end})..."
  chunk_file="${acc_dir}/${taxon}_chunk$(printf '%03d' "$chunk_num").txt"

  cat "$search_xml" \
    | efetch -format acc -start "$start" -stop "$end" \
    > "$chunk_file"

  n_saved=$(grep -c '[^[:space:]]' "$chunk_file" || true)
  echo "    -> ${n_saved} accessions"
  chunk_files+=("$chunk_file")

  start=$(( start + CHUNK_SIZE ))
  chunk_num=$(( chunk_num + 1 ))
  if [[ "$start" -lt "$total" ]]; then
    sleep 1
  fi
done

# ---- concatenate chunks ------------------------------------------------
acc_file="${acc_dir}/${taxon}_ncbi_acc.txt"
cat "${chunk_files[@]}" | grep '[^[:space:]]' > "$acc_file"
rm -f "${chunk_files[@]}"

n_total=$(wc -l < "$acc_file" | tr -d ' ')
echo "  Accession list: ${n_total} entries -> ${acc_file}"

# =======================================================================
# STEP 3 — build the subset BLAST database
# =======================================================================
echo ""
echo "--- Step 3: Building subset BLAST database ---"
echo "  cmd: bash '${subset_script}' '${db_path}' '${acc_file}' '${out_dir}' '${title}'"
bash "$subset_script" "$db_path" "$acc_file" "$out_dir" "$title"

echo ""
echo "=== Done ==="
echo "  Accession list: ${acc_file}"
echo "  Database:       ${out_dir}/$(basename "$out_dir").*"
