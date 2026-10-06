#!/usr/bin/env bash
# make_subset_blastdb.sh
# Fetch accession numbers from NCBI nuccore and build a standalone subset BLAST database.
# Bash equivalent of make_subset_blastdb.R — requires only BLAST+ and curl (no R, no edirect).
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

# ---- load the user's shell PATH from ~/.bashrc ---------------------------
# RStudio's Terminal does not inherit the PATH set up for the macOS Terminal app,
# so BLAST+ may not be found. Load ~/.bashrc before anything else
# (and before strict mode below, since rc files often use unset variables).
bashrc_loaded=0
if [[ -f "$HOME/.bashrc" ]]; then
  _start_dir="$PWD"
  source "$HOME/.bashrc" >/dev/null 2>&1 || true
  cd "$_start_dir"
  bashrc_loaded=1
fi

set -euo pipefail

# record the exact top-level command for provenance (saved in <output_dir>/code/)
wrapper_path="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
wrapper_cmd="$(printf '%q ' bash "$wrapper_path" "$@")"

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

# ---- setup check (printed on every run, as a reminder for new users) -------
print_setup_check() {
  local blast; blast="$(command -v blastdbcmd || true)"
  echo "------------------------------------------------------------------------"
  echo "Setup check. Settings belong in ~/.bashrc, which these scripts load"
  echo "automatically (RStudio does not inherit your Terminal's settings)."
  if [[ "$bashrc_loaded" -eq 1 ]]; then echo "  ~/.bashrc     : loaded"
  else echo "  ~/.bashrc     : NOT FOUND -- create it and add the lines below"; fi
  if [[ -n "${NCBI_EMAIL:-}" ]]; then echo "  NCBI_EMAIL    : ${NCBI_EMAIL}"
  else echo "  NCBI_EMAIL    : not set -- NCBI asks users to identify themselves; add"
       echo "                  export NCBI_EMAIL=\"you@uw.edu\""; fi
  if [[ -n "${NCBI_API_KEY:-}" ]]; then echo "  NCBI_API_KEY  : set"
  else echo "  NCBI_API_KEY  : not set -- recommended (faster, fewer errors); get a free key"
       echo "                  under NCBI account > Account settings > API Key Management, then add"
       echo "                  export NCBI_API_KEY=\"your_key\""; fi
  if [[ -n "$blast" ]]; then echo "  BLAST+        : ${blast} ($(blastdbcmd -version 2>&1 | head -1 | sed 's/^blastdbcmd: //'))"
  else echo "  BLAST+        : NOT FOUND on PATH -- add its bin folder to PATH in ~/.bashrc"; fi
  echo "------------------------------------------------------------------------"
}
print_setup_check

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

# ---- check the source BLAST database before the (possibly long) download --
db_dir="$(dirname "$db_path")"
if [[ ! -d "$db_dir" ]]; then
  echo "ERROR: Source database folder not found: ${db_dir} (is the drive mounted?)" >&2
  exit 1
fi
if ! db_check=$(cd "$db_dir" && blastdbcmd -db "$(basename "$db_path")" -info 2>&1); then
  echo "ERROR: Cannot open the source BLAST database: ${db_path}" >&2
  printf '  blastdbcmd said: %s\n' "$(printf '%s\n' "$db_check" | head -2)" >&2
  db_diagnostics "$db_dir" "$(basename "$db_path")" >&2
  exit 1
fi
echo "Source database OK: ${db_path}"

# ---- NCBI E-utilities settings -------------------------------------------
# Accessions are downloaded directly from NCBI E-utilities with curl, in chunks,
# with retries. NCBI's servers often return transient errors (e.g. HTTP 502) on
# large downloads; each chunk is retried with increasing waits, every chunk is
# checked for completeness, and an interrupted run resumes where it stopped
# when re-run with the same arguments.
# Optional environment variables (e.g. set in ~/.bashrc):
#   NCBI_EMAIL       your email; NCBI asks E-utilities users to identify themselves
#   NCBI_API_KEY     NCBI API key; allows faster requests (free from your NCBI account)
#   NCBI_CHUNK_SIZE  accessions per request (default 5000; NCBI returns at most 9999)
#   NCBI_MAX_TRIES   attempts per request before giving up (default 8)
#   NCBI_RETRY_WAIT  seconds before the first retry; doubles each time (default 5)
RETRY_WAIT="${NCBI_RETRY_WAIT:-5}"
EUTILS="${EUTILS_BASE:-https://eutils.ncbi.nlm.nih.gov/entrez/eutils}"
CHUNK_SIZE="${NCBI_CHUNK_SIZE:-5000}"
MAX_TRIES="${NCBI_MAX_TRIES:-8}"
api_args=(--data-urlencode "tool=SubsetBlastDatabase")   # identifies this tool to NCBI
if [[ -n "${NCBI_EMAIL:-}" ]]; then api_args+=(--data-urlencode "email=${NCBI_EMAIL}"); fi
pause=0.4                                   # NCBI allows 3 requests/s without a key
if [[ -n "${NCBI_API_KEY:-}" ]]; then
  api_args+=(--data-urlencode "api_key=${NCBI_API_KEY}")
  pause=0.15                                # 10 requests/s with a key
fi

# POST one E-utilities request (no retries here); writes the response to $1.
# usage: eutil_post <out_file> <utility.fcgi> [--data-urlencode name=value ...]
eutil_post() {
  local out="$1" util="$2"; shift 2
  curl -sS -f --max-time 300 -X POST "${EUTILS}/${util}" "$@" \
    ${api_args[@]+"${api_args[@]}"} -o "$out"
}

backoff() {   # wait 5, 10, 20, 40 ... seconds (max 300) before attempt $1+1
  local w=$(( RETRY_WAIT * 2 ** ($1 - 1) )); (( w > 300 )) && w=300
  echo "$w"
}

# Run esearch with history; sets globals: total, query_key, web_env
run_esearch() {
  local xml try reason w
  xml="$(mktemp "${TMPDIR:-/tmp}/esearch_XXXXXX")"
  for (( try = 1; try <= MAX_TRIES; try++ )); do
    if eutil_post "$xml" esearch.fcgi --data-urlencode "db=nuccore" \
         --data-urlencode "usehistory=y" --data-urlencode "retmax=0" \
         --data-urlencode "term=${query}" 2>"${xml}.err"; then
      total=$(grep -o '<Count>[0-9]*</Count>' "$xml" | head -1 | tr -dc '0-9' || true)
      query_key=$(grep -o '<QueryKey>[0-9]*</QueryKey>' "$xml" | head -1 | tr -dc '0-9' || true)
      web_env=$(grep -o '<WebEnv>[^<]*</WebEnv>' "$xml" | head -1 | sed 's:</*WebEnv>::g' || true)
      if [[ -n "$total" && -n "$query_key" && -n "$web_env" ]]; then
        rm -f "$xml" "${xml}.err"; return 0
      fi
      reason="unexpected response: $(head -c 300 "$xml" | tr '\n' ' ')"
    else
      reason="$(tail -1 "${xml}.err")"
    fi
    if (( try < MAX_TRIES )); then
      w=$(backoff "$try")
      echo "    esearch attempt ${try}/${MAX_TRIES} failed (${reason}); retrying in ${w}s" >&2
      sleep "$w"
    fi
  done
  rm -f "$xml" "${xml}.err"
  echo "ERROR: esearch failed after ${MAX_TRIES} attempts (${reason})" >&2
  return 1
}

# =======================================================================
# STEP 1 — esearch
# =======================================================================
echo ""
echo "--- Step 1: Querying NCBI nuccore ---"
echo "  query: ${query}"
run_esearch || exit 1
search_date="$(date '+%Y-%m-%d %H:%M %Z')"
if [[ "$total" -eq 0 ]]; then
  echo "ERROR: No records match the query." >&2
  exit 1
fi
search_total="$total"
echo "  Total matching records: ${total}"

# =======================================================================
# STEP 2 — efetch accessions in chunks (with retries and resume)
# =======================================================================
echo ""
echo "--- Step 2: Fetching accessions ---"

acc_dir="${out_dir}/accessions"
mkdir -p "$acc_dir"

# resume support: keep finished chunks only if the search is unchanged
state_file="${acc_dir}/.download_state"
# each request overlaps the next chunk by a few records, so an off-by-one at a chunk
# boundary can never lose an accession. Duplicates from the overlap are removed when
# the chunks are combined.
OVERLAP=10
# NCBI returns at most 9,999 records per efetch request, however many are asked for,
# so keep each request (chunk + overlap) under that limit
NCBI_MAX_PER_REQUEST=9999
if (( CHUNK_SIZE + OVERLAP > NCBI_MAX_PER_REQUEST )); then
  CHUNK_SIZE=$(( NCBI_MAX_PER_REQUEST - OVERLAP ))
  echo "  Note: chunk size reduced to ${CHUNK_SIZE} (NCBI returns at most ${NCBI_MAX_PER_REQUEST} records per request)"
fi
state="$(printf 'query=%s\ntotal=%s\nchunk_size=%s\noverlap=%s' "$query" "$search_total" "$CHUNK_SIZE" "$OVERLAP")"
if [[ -f "$state_file" && "$(cat "$state_file")" == "$state" ]]; then
  echo "  Resuming an earlier download: completed chunks will be reused"
else
  rm -f "${acc_dir}/${taxon}"_chunk*.txt "${acc_dir}/${taxon}"_chunk*.part
  printf '%s\n' "$state" > "$state_file"
fi

# an accession: letters/digits/underscore, optional .version (e.g. MZ463940.1, NC_012920.1,
# or PDB-derived records such as 9SL3_A); anything else (HTML, error text) is rejected
ACC_RE='^[A-Za-z0-9][A-Za-z0-9_.-]*$'
n_chunks=$(( (search_total + CHUNK_SIZE - 1) / CHUNK_SIZE ))
n_short=0
chunk_files=()
for (( i = 1; i <= n_chunks; i++ )); do
  start=$(( (i - 1) * CHUNK_SIZE ))
  want=$(( search_total - start )); (( want > CHUNK_SIZE )) && want=$CHUNK_SIZE
  chunk_file="${acc_dir}/${taxon}_chunk$(printf '%05d' "$i").txt"
  chunk_files+=("$chunk_file")
  label="Chunk ${i}/${n_chunks} (records $(( start + 1 ))-$(( start + want )))"

  # chunk files are only written once complete, so an existing one can be reused
  if [[ -s "$chunk_file" ]]; then
    echo "  ${label}: already downloaded"
    continue
  fi

  ok=0; prev_n=-1; got=0
  for (( try = 1; try <= MAX_TRIES; try++ )); do
    part="${chunk_file}.part"
    if eutil_post "$part" efetch.fcgi --data-urlencode "db=nuccore" \
         --data-urlencode "query_key=${query_key}" --data-urlencode "WebEnv=${web_env}" \
         --data-urlencode "retstart=${start}" --data-urlencode "retmax=$(( want + OVERLAP ))" \
         --data-urlencode "rettype=acc" --data-urlencode "retmode=text" 2>"${part}.err"; then
      n=$(grep -c '[^[:space:]]' "$part" || true)
      bad=$(grep '[^[:space:]]' "$part" | grep -cvE "$ACC_RE" || true)
      # accept a full chunk (>= want, thanks to the overlap), or a short one that NCBI
      # returns identically twice in a row (e.g. the final chunk, or records counted by
      # esearch but never returned); a truncated transfer gives a different count each time
      if [[ "$bad" -eq 0 ]] && { [[ "$n" -ge "$want" ]] || [[ "$n" -gt 0 && "$n" -eq "$prev_n" ]]; }; then
        grep '[^[:space:]]' "$part" > "$chunk_file"
        rm -f "$part" "${part}.err"
        ok=1; got=$n; break
      fi
      reason="got ${n} of ${want} accessions"
      if [[ "$bad" -gt 0 ]]; then
        reason="${reason}; ${bad} line(s) not accessions: $(grep '[^[:space:]]' "$part" | grep -vE "$ACC_RE" | head -1 | cut -c1-120)"
        prev_n=-1
      else
        prev_n=$n
      fi
    else
      reason="$(tail -1 "${part}.err")"
    fi
    if (( try < MAX_TRIES )); then
      w=$(backoff "$try")
      echo "    ${label}: attempt ${try}/${MAX_TRIES} failed (${reason}); retrying in ${w}s" >&2
      sleep "$w"
      # every 3rd failure, start a fresh NCBI search session in case the old one expired
      if (( try % 3 == 0 )); then
        echo "    refreshing the NCBI search session..." >&2
        if run_esearch && [[ "$total" -ne "$search_total" ]]; then
          echo "    WARNING: NCBI now reports ${total} records (was ${search_total}); the database changed during the download." >&2
        fi
      fi
    fi
  done
  rm -f "${chunk_file}.part" "${chunk_file}.part.err"
  if [[ "$ok" -ne 1 ]]; then
    echo "ERROR: ${label} failed after ${MAX_TRIES} attempts (${reason})." >&2
    echo "  Completed chunks are kept in ${acc_dir}; re-run the same command to resume." >&2
    exit 1
  fi
  if [[ "$got" -gt "$want" ]]; then got=$want; fi   # extra lines are overlap
  if [[ "$got" -lt "$want" ]]; then
    echo "  ${label}: ${got} accessions (NCBI consistently returns $(( want - got )) fewer than its count; accepted)"
    n_short=$(( n_short + want - got ))
  else
    echo "  ${label}: ${got} accessions"
  fi
  sleep "$pause"
done

# ---- concatenate chunks (dropping any duplicates) ----------------------
acc_file="${acc_dir}/${taxon}_ncbi_acc.txt"
cat "${chunk_files[@]}" | awk '!seen[$0]++' > "$acc_file"
rm -f "${chunk_files[@]}" "$state_file"
n_total=$(grep -c . "$acc_file")
echo "  Accession list: ${n_total} unique entries -> ${acc_file}"
if [[ "$n_total" -ne "$search_total" ]]; then
  echo "  Note: NCBI's search counted ${search_total} records; ${n_total} unique accessions were retrieved." >&2
  echo "        (NCBI's count includes some records it does not return, e.g. suppressed or" >&2
  echo "        withdrawn ones; records can also change at NCBI during a long download.)" >&2
fi
total="$search_total"

# =======================================================================
# STEP 3 — build the subset BLAST database
# =======================================================================
echo ""
echo "--- Step 3: Building subset BLAST database ---"
echo "  cmd: bash '${subset_script}' '${db_path}' '${acc_file}' '${out_dir}' '${title}'"
SUBSET_SEARCH_QUERY="$query" \
SUBSET_SEARCH_DATE="$search_date" \
SUBSET_SEARCH_COUNT="$total" \
SUBSET_CALLER_SCRIPT="$wrapper_path" \
SUBSET_CALLER_COMMAND="$wrapper_cmd" \
  bash "$subset_script" "$db_path" "$acc_file" "$out_dir" "$title"

echo ""
echo "=== Done ==="
echo "  Accession list: ${acc_file}"
echo "  readme.txt:     ${out_dir}/readme.txt"
echo "  Database:       ${out_dir}/$(printf '%s' "$(basename "$out_dir")" | tr -c 'A-Za-z0-9._-' '_').*"
