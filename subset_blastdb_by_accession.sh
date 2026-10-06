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
#
# Provenance: a copy of this script, the exact command used, and the input
# accession list are saved to <output_dir>/code/, and a short README.txt
# describing the database contents and dates is written to <output_dir>/.
# When called from make_subset_blastdb.sh, the NCBI search query, date and
# record count (passed via SUBSET_SEARCH_* env vars) are recorded as well,
# and the wrapper script is copied into code/ too.

# ---- load the user's shell PATH from ~/.bashrc ---------------------------
# RStudio's Terminal does not inherit the PATH set up for the macOS Terminal app,
# so BLAST+ and edirect may not be found. Load ~/.bashrc before anything else
# (and before strict mode below, since rc files often use unset variables).
bashrc_loaded=0
if [[ -f "$HOME/.bashrc" ]]; then
  _start_dir="$PWD"
  source "$HOME/.bashrc" >/dev/null 2>&1 || true
  cd "$_start_dir"
  bashrc_loaded=1
fi

set -euo pipefail

# ---- record provenance before args are consumed ------------------------
script_path="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
invocation="$(printf '%q ' bash "$script_path" "$@")"
run_dir="$(pwd)"

# ---- parse args --------------------------------------------------------
if [[ $# -lt 2 ]]; then
  echo "Usage: bash subset_blastdb_by_accession.sh <blastdb> <acc_list> [output_dir] [title]" >&2
  exit 1
fi

db_path="$1"
acc_file="$2"
out_dir="${3:-subset_db}"
prefix="$(basename "$out_dir")"

mkdir -p "$out_dir"
out_abs="$(cd "$out_dir" && pwd -P)"   # absolute path; used for all intermediate files

# ---- BLAST+ and paths with spaces --------------------------------------
# BLAST+ splits database names on spaces (-db "a b" means two databases), so a
# database path under e.g. "Kelly_Lab - Documents/2. KellyLab" breaks blastdbcmd
# and makeblastdb. Workaround: run each BLAST command from inside the database's
# own folder and refer to the database by its bare (space-free) name.
if [[ ! -d "$(dirname "$db_path")" ]]; then
  echo "ERROR: Source database folder not found: $(dirname "$db_path") (is the drive mounted?)" >&2
  exit 1
fi
db_dir="$(cd "$(dirname "$db_path")" && pwd -P)"
db_name="$(basename "$db_path")"
if [[ "$db_name" =~ [[:space:]] ]]; then
  echo "ERROR: source database name contains spaces, which BLAST+ cannot handle: ${db_name}" >&2
  exit 1
fi

# file prefix from the output dir name; characters BLAST+ can't handle become "_"
prefix="$(printf '%s' "$prefix" | tr -c 'A-Za-z0-9._-' '_')"
db_title="${4:-$prefix}"
fasta_file="${out_abs}/${prefix}.fasta"
taxid_file="${out_abs}/${prefix}_taxids.txt"
taxid_map="${out_abs}/${prefix}_taxid_map.txt"
new_db_name="${out_abs}/${prefix}"

# ---- 0. validate input -------------------------------------------------
# Count non-blank lines
n_accs=$(grep -c '[^[:space:]]' "$acc_file" || true)
echo "Input: ${n_accs} accession(s) from ${acc_file}"

# Warn if any lines look like raw GI numbers (pure integers)
gi_count=$(grep -cE '^[0-9]+$' "$acc_file" 2>/dev/null || true)
if [[ "$gi_count" -gt 0 ]]; then
  echo "WARNING: ${gi_count} line(s) look like GI numbers (pure integers). Verify your accession list." >&2
fi

# Explain why a BLAST database could not be opened (shown with the error)
db_diagnostics() {   # $1 = db folder, $2 = db name
  local d="$1" n="$2" listing idx
  echo "  Diagnostics:"
  echo "    blastdbcmd used : $(command -v blastdbcmd || echo 'NOT FOUND') ($(blastdbcmd -version 2>&1 | head -1))"
  if ! listing=$(ls "$d" 2>&1); then
    echo "    cannot list ${d}: $(printf '%s\n' "$listing" | head -1)"
    echo "    If that says 'Operation not permitted', macOS is blocking this app from the drive:"
    echo "    allow RStudio (and/or Terminal) under System Settings > Privacy & Security >"
    echo "    Files and Folders (Removable Volumes), or Full Disk Access, then restart it."
  else
    idx=$(printf '%s\n' "$listing" | grep -E "^${n//./\\.}(\.[0-9]+)?\.(nal|nin|ndb)$" | head -4 | tr '\n' ' ' || true)
    if [[ -n "$idx" ]]; then
      echo "    '${n}' index files in ${d}: ${idx}"
      echo "    The files are there, so the BLAST+ version above may be too old for this database."
    else
      echo "    '${n}' index files in ${d}: NONE"
      echo "    The folder holds $(printf '%s\n' "$listing" | grep -c .) items, e.g.: $(printf '%s\n' "$listing" | head -4 | tr '\n' ' ')"
      echo "    Check the database name in the path."
    fi
  fi
}

# ---- check the source database can be opened before doing any work ----
if ! db_info=$(cd "$db_dir" && blastdbcmd -db "$db_name" -info 2>&1); then
  rmdir "$out_dir" 2>/dev/null || true   # remove output dir only if we just made it empty
  {
    echo "ERROR: Cannot open the source BLAST database: ${db_path}"
    echo "  blastdbcmd said: $(printf '%s\n' "$db_info" | head -3)"
    db_diagnostics "$db_dir" "$db_name"
    echo "  (An LMDB / mdb_env_open error usually means the drive is network-mounted or the"
    echo "   database download is incomplete.)"
  } >&2
  exit 1
fi
echo "Source database opened OK: ${db_path}"

# ---- detect available cores --------------------------------------------
if command -v nproc &>/dev/null; then
  n_cores=$(nproc)
else
  n_cores=$(sysctl -n hw.logicalcpu 2>/dev/null || echo 1)
fi
echo "Parallel workers: ${n_cores}"

# ---- split accession list into N chunks --------------------------------
chunk_prefix="${out_abs}/.${prefix}_acc_chunk_"
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
  ( cd "$db_dir" && blastdbcmd \
      -db "$db_name" \
      -entry_batch "$chunk" \
      -outfmt $'>%a %t\n%s' \
      -out "$cf" 2> "${cf}.err" || true ) &
done
wait
cat "${fasta_chunks[@]}" > "$fasta_file"
# blastdbcmd reports each accession missing from the source db as "Skipped <acc>";
# collect those in a file instead of printing thousands of lines
not_found_file="${out_abs}/${prefix}_not_in_source_db.txt"
cat "${fasta_chunks[@]/%/.err}" | sed -n 's/.*Skipped[[:space:]]*//p' > "$not_found_file"
n_not_found=$(grep -c . "$not_found_file" || true)
other_errs=$(cat "${fasta_chunks[@]/%/.err}" | grep -v 'Skipped' | grep -v 'Entry or entries not found' | grep . | head -5 || true)
[[ -n "$other_errs" ]] && printf '  blastdbcmd messages:\n%s\n' "$other_errs" >&2
rm -f "${fasta_chunks[@]}" "${fasta_chunks[@]/%/.err}"
if [[ "$n_not_found" -gt 0 ]]; then
  echo "  ${n_not_found} accession(s) not in the source database (listed in $(basename "$not_found_file"))"
else
  rm -f "$not_found_file"
fi
if [[ ! -s "$fasta_file" ]]; then
  echo "ERROR: blastdbcmd (fasta) produced no output: none of the accessions were found in ${db_path}" >&2
  rm -f "${acc_chunks[@]}" "$fasta_file"
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
  ( cd "$db_dir" && blastdbcmd \
      -db "$db_name" \
      -entry_batch "$chunk" \
      -outfmt "%a	%T" \
      -out "$ct" 2>/dev/null || true ) &
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
echo "  cmd: (cd '${out_abs}' && makeblastdb -in '${prefix}.fasta' -dbtype nucl -parse_seqids -taxid_map '${prefix}_taxid_map.txt' -out '${prefix}' -title '${db_title}')"
# run inside out_abs with bare file names (see "BLAST+ and paths with spaces" above)
( cd "$out_abs" && makeblastdb \
    -in "$(basename "$fasta_file")" \
    -dbtype nucl \
    -parse_seqids \
    -taxid_map "$(basename "$taxid_map")" \
    -out "$prefix" \
    -title "$db_title" )

# ---- 6. save code + write readme --------------------------------------
n_taxids=$(awk '{print $2}' "$taxid_map" | sort -u | wc -l | tr -d ' ')

# copy of the generating code, the exact command, and the input accession list
code_dir="${out_dir}/code"
mkdir -p "$code_dir"
cp "$script_path" "$code_dir/"
code_line="code/$(basename "$script_path")"
caller_script="${SUBSET_CALLER_SCRIPT:-}"
if [[ -n "$caller_script" && -f "$caller_script" ]]; then
  cp "$caller_script" "$code_dir/"
  code_line="code/$(basename "$caller_script") (calls ${code_line})"
fi

# keep the accession list: reference it if already inside out_dir, else copy it
acc_abs="$(cd "$(dirname "$acc_file")" && pwd -P)/$(basename "$acc_file")"
out_abs="$(cd "$out_dir" && pwd -P)"
if [[ "$acc_abs" == "${out_abs}/"* ]]; then
  acc_rel="${acc_abs#"${out_abs}/"}"
else
  cp "$acc_file" "${code_dir}/input_accessions.txt"
  acc_rel="code/input_accessions.txt (copy of $(basename "$acc_file"))"
fi

{
  echo "# Run from: ${run_dir}"
  echo "# Run on:   $(date '+%Y-%m-%d %H:%M:%S %Z') by $(whoami)@$(hostname)"
  echo "# BLAST+:   $(blastdbcmd -version 2>/dev/null | head -1 || echo unknown)"
  if [[ -n "${SUBSET_CALLER_COMMAND:-}" ]]; then
    echo "# Top-level command:"
    echo "${SUBSET_CALLER_COMMAND}"
    echo "# which ran:"
  fi
  echo "${invocation}"
} > "${code_dir}/command.txt"

# date of the source database snapshot (from blastdbcmd -info)
src_db_date=$(printf '%s\n' "$db_info" \
  | sed -n 's/.*Date: \(.*[AP]M\).*/\1/p' | head -1 || true)
[[ -z "$src_db_date" ]] && src_db_date="unknown"

readme_file="${out_dir}/readme.txt"
{
  echo "${db_title}"
  echo "Standalone nucleotide BLAST database (${prefix}.*), built $(date '+%Y-%m-%d') by $(whoami)."
  echo "Contents: ${n_after} unique accessions, ${n_taxids} unique NCBI taxon IDs,"
  echo "  from the accession list ${acc_rel}."
  if [[ "$n_not_found" -gt 0 ]]; then
    echo "  ${n_not_found} listed accessions were not in the source database: $(basename "$not_found_file")."
  fi
  if [[ -n "${SUBSET_SEARCH_QUERY:-}" ]]; then
    echo "Search: NCBI nuccore, run ${SUBSET_SEARCH_DATE:-unknown date}, ${SUBSET_SEARCH_COUNT:-?} records matched:"
    echo "  ${SUBSET_SEARCH_QUERY}"
  fi
  echo "Source: ${db_path} (source db dated ${src_db_date})."
  echo "Code: ${code_line}; exact command(s) in code/command.txt."
} > "$readme_file"
echo "readme.txt written: ${readme_file}"

# ---- 7. clean up -------------------------------------------------------
rm "$fasta_file"
echo "Removed intermediate FASTA: ${fasta_file}"

echo ""
echo "Done. New database: ${new_db_name}"
echo "Files created:"
echo "  readme.txt: ${readme_file}"
echo "  Code:       ${code_dir}/"
echo "  Taxid list: ${taxid_file}"
echo "  Taxid map:  ${taxid_map}"
echo "  BLAST db:   ${new_db_name}.*"
