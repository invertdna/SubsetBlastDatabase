#!/usr/bin/env bash
# subset_blastdb_by_accession.sh
# Create a standalone nucleotide BLAST database from a subset of NCBI accession numbers.
# Bash equivalent of subset_blastdb_by_accession.R — requires only BLAST+ (no R).
#
# Usage:
#   bash subset_blastdb_by_accession.sh <path_to_blastdb> <path_to_acc_list> [output_dir] [title]
#
# Arguments:
#   path_to_blastdb : path (with db name prefix) to an existing local BLAST db
#                     e.g. /Volumes/Clupea/core_nt/core_nt
#   path_to_acc_list: text file with one accession (or accession.version) per line
#   output_dir      : (optional) directory for output files; default "subset_db"
#   title           : (optional) human-readable title for the database; default = output_dir name
#
# The accession list is split into N chunks and blastdbcmd is run in parallel,
# where N = number of logical CPU cores. Chunk files are removed after merging.

set -euo pipefail

# ---- parse args --------------------------------------------------------
if [[ $# -lt 2 ]]; then
  echo "Usage: bash subset_blastdb_by_accession.sh <blastdb> <acc_list> [output_dir] [title]" >&2
  exit 1
fi

db_path="$1"
acc_file="$2"
out_dir="${3:-subset_db}"
prefix="$(basename "$out_dir")"
db_title="${4:-$prefix}"

mkdir -p "$out_dir"

fasta_file="${out_dir}/${prefix}.fasta"
taxid_file="${out_dir}/${prefix}_taxids.txt"
taxid_map="${out_dir}/${prefix}_taxid_map.txt"
new_db_name="${out_dir}/${prefix}"

# ---- 0. validate input -------------------------------------------------
# Count non-blank lines
n_accs=$(grep -c '[^[:space:]]' "$acc_file" || true)
echo "Input: ${n_accs} accession(s) from ${acc_file}"

# Warn if any lines look like raw GI numbers (pure integers)
gi_count=$(grep -cE '^[0-9]+$' "$acc_file" 2>/dev/null || true)
if [[ "$gi_count" -gt 0 ]]; then
  echo "WARNING: ${gi_count} line(s) look like GI numbers (pure integers). Verify your accession list." >&2
fi

# ---- detect available cores --------------------------------------------
if command -v nproc &>/dev/null; then
  n_cores=$(nproc)
else
  n_cores=$(sysctl -n hw.logicalcpu 2>/dev/null || echo 1)
fi
echo "Parallel workers: ${n_cores}"

# ---- split accession list into N chunks --------------------------------
chunk_prefix="${out_dir}/.${prefix}_acc_chunk_"
lines_per_chunk=$(( (n_accs + n_cores - 1) / n_cores ))
[[ "$lines_per_chunk" -lt 1 ]] && lines_per_chunk=1
split -l "$lines_per_chunk" -d "$acc_file" "$chunk_prefix"
acc_chunks=( "${chunk_prefix}"* )
n_chunks=${#acc_chunks[@]}
echo "Split into ${n_chunks} chunk(s) of up to ${lines_per_chunk} accessions"

# ---- 1. extract FASTA in parallel --------------------------------------
echo "Extracting FASTA sequences (${n_chunks} parallel job(s))..."
fasta_chunks=()
for chunk in "${acc_chunks[@]}"; do
  cf="${chunk}.fasta"
  : > "$cf"   # ensure file exists even if blastdbcmd finds nothing
  fasta_chunks+=("$cf")
  ( blastdbcmd \
      -db "$db_path" \
      -entry_batch "$chunk" \
      -outfmt $'>%a %t\n%s' \
      -out "$cf" || true ) &
done
wait
cat "${fasta_chunks[@]}" > "$fasta_file"
rm "${fasta_chunks[@]}"
if [[ ! -s "$fasta_file" ]]; then
  echo "ERROR: blastdbcmd (fasta) produced no output" >&2
  rm "${acc_chunks[@]}" 2>/dev/null || true
  exit 1
fi
echo "  FASTA chunks merged"

# ---- 2. extract taxids in parallel -------------------------------------
echo "Extracting taxon IDs (${n_chunks} parallel job(s))..."
taxid_chunks=()
for chunk in "${acc_chunks[@]}"; do
  ct="${chunk}.taxids"
  : > "$ct"   # ensure file exists even if blastdbcmd finds nothing
  taxid_chunks+=("$ct")
  ( blastdbcmd \
      -db "$db_path" \
      -entry_batch "$chunk" \
      -outfmt "%a	%T" \
      -out "$ct" || true ) &
done
wait
cat "${taxid_chunks[@]}" > "$taxid_file"
rm "${taxid_chunks[@]}" "${acc_chunks[@]}"
if [[ ! -s "$taxid_file" ]]; then
  echo "ERROR: blastdbcmd (taxid) produced no output" >&2
  exit 1
fi
echo "  Taxid chunks merged"

# ---- 3. build deduplicated taxid map -----------------------------------
# Keep first occurrence per accession; drop taxid == 0 or empty
awk -F'\t' '!seen[$1]++ && $2 != "" && $2 != "0"' "$taxid_file" > "$taxid_map"
n_map=$(wc -l < "$taxid_map" | tr -d ' ')
echo "Taxid map written: ${n_map} entries"

# ---- 4. deduplicate FASTA ----------------------------------------------
# blastdbcmd -outfmt "%s" writes each sequence as a single line, so every
# record is exactly 2 lines: a > header and a sequence line.
# Keep only the first occurrence of each accession.
n_before=$(grep -c '^>' "$fasta_file" || true)
awk '
  /^>/ {
    acc = substr($1, 2)
    if (acc in seen) { skip = 1 } else { seen[acc] = 1; skip = 0 }
  }
  !skip { print }
' "$fasta_file" > "${fasta_file}.tmp" && mv "${fasta_file}.tmp" "$fasta_file"
n_after=$(grep -c '^>' "$fasta_file" || true)
n_dedup=$(( n_before - n_after ))
if [[ "$n_dedup" -gt 0 ]]; then
  echo "Removed ${n_dedup} duplicate accession(s) from FASTA"
fi

# ---- 5. build BLAST database -------------------------------------------
echo "Building new BLAST database (title=${db_title})..."
echo "  cmd: makeblastdb -in '${fasta_file}' -dbtype nucl -parse_seqids -taxid_map '${taxid_map}' -out '${new_db_name}' -title '${db_title}'"
makeblastdb \
  -in "$fasta_file" \
  -dbtype nucl \
  -parse_seqids \
  -taxid_map "$taxid_map" \
  -out "$new_db_name" \
  -title "$db_title"

# ---- 6. write readme ---------------------------------------------------
n_taxids=$(awk '{print $2}' "$taxid_map" | sort -u | wc -l | tr -d ' ')
readme_file="${out_dir}/readme.txt"
cat > "$readme_file" <<EOF
SubsetBlastDatabase
===================
Script:            $(basename "$0")
Source database:   ${db_path}
Created by:        $(whoami)
Created:           $(date)
Unique accessions: ${n_after}
Unique taxon IDs:  ${n_taxids}
EOF
echo "readme.txt written: ${readme_file}"

# ---- 7. clean up -------------------------------------------------------
rm "$fasta_file"
echo "Removed intermediate FASTA: ${fasta_file}"

echo ""
echo "Done. New database: ${new_db_name}"
echo "Files created:"
echo "  readme.txt: ${readme_file}"
echo "  Taxid list: ${taxid_file}"
echo "  Taxid map:  ${taxid_map}"
echo "  BLAST db:   ${new_db_name}.*"
