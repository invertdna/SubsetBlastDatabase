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

# ---- 1. extract FASTA --------------------------------------------------
echo "Extracting FASTA sequences..."
echo "  cmd: blastdbcmd -db '${db_path}' -entry_batch '${acc_file}' -outfmt '>%a %t\n%s' -out '${fasta_file}'"
blastdbcmd \
  -db "$db_path" \
  -entry_batch "$acc_file" \
  -outfmt $'>%a %t\n%s' \
  -out "$fasta_file" || {
    status=$?
    if [[ ! -s "$fasta_file" ]]; then
      echo "ERROR: blastdbcmd (fasta) produced no output (exit code ${status})" >&2
      exit "$status"
    fi
    echo "  Note: blastdbcmd skipped some accessions not present in local db (exit code ${status})"
  }

# ---- 2. extract taxids -------------------------------------------------
echo "Extracting taxon IDs..."
echo "  cmd: blastdbcmd -db '${db_path}' -entry_batch '${acc_file}' -outfmt '%a\t%T' -out '${taxid_file}'"
blastdbcmd \
  -db "$db_path" \
  -entry_batch "$acc_file" \
  -outfmt "%a	%T" \
  -out "$taxid_file" || {
    status=$?
    if [[ ! -s "$taxid_file" ]]; then
      echo "ERROR: blastdbcmd (taxid) produced no output (exit code ${status})" >&2
      exit "$status"
    fi
    echo "  Note: blastdbcmd skipped some accessions not present in local db (exit code ${status})"
  }

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

# ---- 6. clean up -------------------------------------------------------
rm "$fasta_file"
echo "Removed intermediate FASTA: ${fasta_file}"

echo ""
echo "Done. New database: ${new_db_name}"
echo "Files created:"
echo "  Taxid list: ${taxid_file}"
echo "  Taxid map:  ${taxid_map}"
echo "  BLAST db:   ${new_db_name}.*"
