#!/usr/bin/env Rscript
# make_subset_blastdb.R
# Fetch accession numbers from NCBI nuccore and build a standalone subset BLAST database.
#
# Usage:
#   Rscript make_subset_blastdb.R <blastdb> <taxon> <query> <title> [output_dir]
#
# Arguments:
#   blastdb    : path (with db prefix) to an existing local BLAST db
#                  e.g. /Volumes/Clupea/blastdb/core_nt
#   taxon      : short name used for intermediate file naming (no spaces)
#                  e.g. Sebastes
#   query      : NCBI esearch query string (quote in shell if it contains spaces)
#                  e.g. "(Sebastes[Organism]) AND (mitochondrion OR mitochondrial) AND 100:20000[Sequence Length]"
#   title      : human-readable title embedded in the output database
#                  e.g. "Sebastes mitochondrial DNA"
#   output_dir : (optional) directory for all output; default = taxon value

# ---- load the user's shell settings from ~/.bashrc ----
# RStudio (including its Terminal and R's system() calls) does not inherit the settings
# made for the macOS Terminal app, so BLAST+ may not be found and NCBI settings are
# missing. Source ~/.bashrc in bash and adopt its PATH, BLASTDB and NCBI_* variables
# for everything this script runs (variables already set in R are kept).
bashrc_loaded <- FALSE
if (file.exists(path.expand("~/.bashrc"))) {
  bashrc_vars <- c("PATH", "BLASTDB", "NCBI_EMAIL", "NCBI_API_KEY",
                   "NCBI_CHUNK_SIZE", "NCBI_MAX_TRIES", "NCBI_RETRY_WAIT")
  bashrc_env <- suppressWarnings(system(paste0(
    "bash -c 'source ~/.bashrc >/dev/null 2>&1; for v in ", paste(bashrc_vars, collapse = " "),
    "; do printf \"%s=%s\\n\" \"$v\" \"${!v}\"; done'"), intern = TRUE))
  for (kv in bashrc_env) {
    k <- sub("=.*", "", kv); v <- sub("^[^=]*=", "", kv)
    if (nzchar(v) && (k == "PATH" || !nzchar(Sys.getenv(k)))) do.call(Sys.setenv, setNames(list(v), k))
  }
  bashrc_loaded <- length(bashrc_env) > 0
}

# ---- parse args -------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop(paste(
    "Usage: Rscript make_subset_blastdb.R <blastdb> <taxon> <query> <title> [output_dir]",
    "  blastdb    path to existing local BLAST db",
    "  taxon      short name for file naming (no spaces)",
    "  query      NCBI esearch query string",
    "  title      title for the output database",
    "  output_dir (optional) output directory; default = taxon",
    sep = "\n"
  ))
}

db_path <- args[1]
taxon   <- args[2]
query   <- args[3]
title   <- args[4]
out_dir <- if (length(args) >= 5) args[5] else taxon

# ---- locate subset_blastdb_by_accession.R next to this script ---------
# ---- locate this script (robust to spaces/special characters in paths) ----
# Rscript encodes spaces in the --file= argument as "~+~", so decode before use.
# Falls back to source()'s ofile, then the RStudio editor path, when not run via Rscript.
this_script_path <- function() {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f) > 0) {
    p <- gsub("~+~", " ", sub("^--file=", "", f[1]), fixed = TRUE)
    return(normalizePath(p, mustWork = TRUE))
  }
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(normalizePath(of, mustWork = TRUE))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- rstudioapi::getSourceEditorContext()$path
    if (nzchar(p)) return(normalizePath(p, mustWork = TRUE))
  }
  NA_character_
}

script_path  <- this_script_path()
script_dir   <- if (!is.na(script_path)) dirname(script_path) else getwd()
subset_script <- file.path(script_dir, "subset_blastdb_by_accession.R")
if (!file.exists(subset_script)) {
  stop("Cannot find subset_blastdb_by_accession.R in: ", script_dir)
}

# ---- helper: run a shell command, stop on failure ---------------------
run <- function(cmd, label = cmd) {
  message("  cmd: ", cmd)
  status <- system(cmd)
  if (status != 0) stop(label, " failed with exit code ", status)
  invisible(status)
}

# ---- setup check (printed on every run, as a reminder for new users) ----
local({
  blast <- Sys.which("blastdbcmd")
  msg <- c(strrep("-", 72),
    "Setup check. Settings belong in ~/.bashrc, which these scripts load",
    "automatically (RStudio does not inherit your Terminal's settings).",
    if (bashrc_loaded) "  ~/.bashrc     : loaded"
    else "  ~/.bashrc     : NOT FOUND -- create it and add the lines below",
    if (nzchar(Sys.getenv("NCBI_EMAIL"))) paste0("  NCBI_EMAIL    : ", Sys.getenv("NCBI_EMAIL"))
    else c("  NCBI_EMAIL    : not set -- NCBI asks users to identify themselves; add",
           '                  export NCBI_EMAIL="you@uw.edu"'),
    if (nzchar(Sys.getenv("NCBI_API_KEY"))) "  NCBI_API_KEY  : set"
    else c("  NCBI_API_KEY  : not set -- recommended (faster, fewer errors); get a free key",
           "                  under NCBI account > Account settings > API Key Management, then add",
           '                  export NCBI_API_KEY="your_key"'),
    if (nzchar(blast)) paste0("  BLAST+        : ", blast, " (",
      sub("^blastdbcmd: ", "", suppressWarnings(system("blastdbcmd -version 2>&1", intern = TRUE))[1]), ")")
    else "  BLAST+        : NOT FOUND on PATH -- add its bin folder to PATH in ~/.bashrc",
    strrep("-", 72))
  message(paste(msg, collapse = "\n"))
})

# Explain why a BLAST database could not be opened (shown with the error)
db_diagnostics <- function(db_dir, db_name) {
  bin <- Sys.which("blastdbcmd")
  ver <- tryCatch(suppressWarnings(system("blastdbcmd -version 2>&1", intern = TRUE))[1],
                  error = function(e) "?")
  out <- c("  Diagnostics:",
           sprintf("    blastdbcmd used : %s (%s)", if (nzchar(bin)) bin else "NOT FOUND", ver))
  ls_out <- suppressWarnings(system(paste("ls", shQuote(db_dir), "2>&1"), intern = TRUE))
  if (!is.null(attr(ls_out, "status")) && attr(ls_out, "status") != 0) {
    out <- c(out, paste0("    cannot list ", db_dir, ": ", ls_out[1]),
      "    If that says 'Operation not permitted', macOS is blocking this app from the drive:",
      "    allow RStudio (and/or Terminal) under System Settings > Privacy & Security >",
      "    Files and Folders (Removable Volumes), or Full Disk Access, then restart RStudio.")
  } else {
    idx <- grep(paste0("^", gsub(".", "\\.", db_name, fixed = TRUE),
                       "(\\.[0-9]+)?\\.(nal|nin|ndb)$"), ls_out, value = TRUE)
    out <- c(out, sprintf("    '%s' index files in %s: %s", db_name, db_dir,
                          if (length(idx)) paste(head(idx, 4), collapse = " ") else "NONE"),
      if (!length(idx)) c(sprintf("    The folder holds %d items, e.g.: %s", length(ls_out),
                                  paste(head(ls_out, 4), collapse = " ")),
                          "    Check the database name in the path.")
      else "    The files are there, so the BLAST+ version above may be too old for this database.")
  }
  paste(out, collapse = "\n")
}

# ---- check the source BLAST database before the (possibly long) download --
if (!dir.exists(dirname(db_path)))
  stop("Source database folder not found: ", dirname(db_path), " (is the drive mounted?)", call. = FALSE)
db_check <- suppressWarnings(system(paste(
  "cd", shQuote(dirname(db_path)), "&& blastdbcmd -db", shQuote(basename(db_path)), "-info 2>&1"),
  intern = TRUE))
if (!is.null(attr(db_check, "status")) && attr(db_check, "status") != 0)
  stop("Cannot open the source BLAST database: ", db_path, "\n  blastdbcmd said: ",
       paste(head(db_check, 2), collapse = "\n    "), "\n",
       db_diagnostics(dirname(db_path), basename(db_path)), call. = FALSE)
message("Source database OK: ", db_path)

# ---- NCBI E-utilities settings ----------------------------------------
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
retry_wait  <- as.numeric(Sys.getenv("NCBI_RETRY_WAIT", "5"))
eutils_base <- Sys.getenv("EUTILS_BASE", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils")
chunk_size  <- as.numeric(Sys.getenv("NCBI_CHUNK_SIZE", "5000"))
max_tries   <- as.integer(Sys.getenv("NCBI_MAX_TRIES", "8"))
api_key     <- Sys.getenv("NCBI_API_KEY")
pause       <- if (nzchar(api_key)) 0.15 else 0.4   # NCBI: 10 req/s with key, 3 without
# an accession: letters/digits/underscore, optional .version (e.g. MZ463940.1, NC_012920.1,
# or PDB-derived records such as 9SL3_A); anything else (HTML, error text) is rejected
acc_regex   <- "^[A-Za-z0-9][A-Za-z0-9_.-]*$"

# POST one E-utilities request (no retries here); returns list(ok, msg)
eutil_post <- function(util, params, out) {
  params$tool <- "SubsetBlastDatabase"                       # identifies this tool to NCBI
  if (nzchar(Sys.getenv("NCBI_EMAIL"))) params$email <- Sys.getenv("NCBI_EMAIL")
  if (nzchar(api_key)) params$api_key <- api_key
  err  <- tempfile()
  args <- c("-sS", "-f", "--max-time", "300", "-X", "POST",
            shQuote(paste0(eutils_base, "/", util)), "-o", shQuote(out))
  for (nm in names(params))
    args <- c(args, "--data-urlencode", shQuote(paste0(nm, "=", params[[nm]])))
  status <- suppressWarnings(system2("curl", args, stdout = FALSE, stderr = err))
  msg <- if (file.exists(err)) tail(c("", readLines(err, warn = FALSE)), 1) else ""
  unlink(err)
  list(ok = identical(as.integer(status), 0L), msg = msg)
}

backoff <- function(try) min(300, retry_wait * 2^(try - 1))   # 5, 10, 20, 40 ... s

# esearch with history; returns list(total, query_key, web_env)
run_esearch <- function() {
  reason <- ""
  for (try in seq_len(max_tries)) {
    out <- tempfile(fileext = ".xml")
    r <- eutil_post("esearch.fcgi",
                    list(db = "nuccore", usehistory = "y", retmax = "0", term = query), out)
    if (r$ok) {
      x <- paste(readLines(out, warn = FALSE), collapse = "\n")
      tag <- function(t) {   # first occurrence (the first <Count> is the overall total)
        m <- regmatches(x, regexpr(sprintf("<%s>[^<]*</%s>", t, t), x))
        if (length(m)) sub("^<[^>]+>([^<]*)</[^>]+>$", "\\1", m) else NA_character_
      }
      res <- list(total = as.numeric(tag("Count")), query_key = tag("QueryKey"),
                  web_env = tag("WebEnv"))
      unlink(out)
      if (!anyNA(unlist(res))) return(res)
      reason <- paste("unexpected response:", substr(gsub("\n", " ", x), 1, 300))
    } else {
      unlink(out)
      reason <- r$msg
    }
    if (try < max_tries) {
      message(sprintf("    esearch attempt %d/%d failed (%s); retrying in %ds",
                      try, max_tries, reason, backoff(try)))
      Sys.sleep(backoff(try))
    }
  }
  stop(sprintf("esearch failed after %d attempts (%s)", max_tries, reason), call. = FALSE)
}

# =======================================================================
# STEP 1 — esearch: get total count and a history session for paging
# =======================================================================
message("\n--- Step 1: Querying NCBI nuccore ---")
message("  query: ", query)
srch <- run_esearch()
search_date <- format(Sys.time(), "%Y-%m-%d %H:%M %Z")
total <- srch$total
message(sprintf("  Total matching records: %.0f", total))
if (total == 0) stop("No records match the query.")

# =======================================================================
# STEP 2 — efetch accessions in chunks (with retries and resume)
# =======================================================================
message("\n--- Step 2: Fetching accessions ---")

acc_dir <- file.path(out_dir, "accessions")
dir.create(acc_dir, showWarnings = FALSE, recursive = TRUE)

# resume support: keep finished chunks only if the search is unchanged
state_file <- file.path(acc_dir, ".download_state")
# each request overlaps the next chunk by a few records, so an off-by-one at a chunk
# boundary can never lose an accession. Duplicates from the overlap are removed when
# the chunks are combined.
overlap <- 10
# NCBI returns at most 9,999 records per efetch request, however many are asked for,
# so keep each request (chunk + overlap) under that limit
ncbi_max_per_request <- 9999
if (chunk_size + overlap > ncbi_max_per_request) {
  chunk_size <- ncbi_max_per_request - overlap
  message(sprintf("  Note: chunk size reduced to %.0f (NCBI returns at most %d records per request)",
                  chunk_size, ncbi_max_per_request))
}
state <- c(paste0("query=", query), sprintf("total=%.0f", total),
           sprintf("chunk_size=%.0f", chunk_size), sprintf("overlap=%d", overlap))
if (file.exists(state_file) && identical(readLines(state_file, warn = FALSE), state)) {
  message("  Resuming an earlier download: completed chunks will be reused")
} else {
  old <- list.files(acc_dir, pattern = paste0("^", gsub("([][{}()+*^$.|\\\\?])", "\\\\\\1", taxon),
                                              "_chunk.*\\.(txt|part)$"), full.names = TRUE)
  invisible(file.remove(old))
  writeLines(state, state_file)
}

n_chunks    <- ceiling(total / chunk_size)
chunk_files <- character(n_chunks)
for (i in seq_len(n_chunks)) {
  start <- (i - 1) * chunk_size
  want  <- min(chunk_size, total - start)
  chunk_file <- file.path(acc_dir, sprintf("%s_chunk%05d.txt", taxon, i))
  chunk_files[i] <- chunk_file
  label <- sprintf("Chunk %d/%d (records %.0f-%.0f)", i, n_chunks, start + 1, start + want)

  # chunk files are only written once complete, so an existing one can be reused
  if (file.exists(chunk_file) && file.size(chunk_file) > 0) {
    message("  ", label, ": already downloaded")
    next
  }

  ok <- FALSE; prev_n <- -1; got <- 0
  for (try in seq_len(max_tries)) {
    part <- paste0(chunk_file, ".part")
    r <- eutil_post("efetch.fcgi",
                    list(db = "nuccore", query_key = srch$query_key, WebEnv = srch$web_env,
                         retstart = sprintf("%.0f", start), retmax = sprintf("%.0f", want + overlap),
                         rettype = "acc", retmode = "text"), part)
    if (r$ok) {
      lines <- trimws(readLines(part, warn = FALSE))
      lines <- lines[nzchar(lines)]
      bad   <- lines[!grepl(acc_regex, lines)]
      n <- length(lines)
      # accept a full chunk (>= want, thanks to the overlap), or a short one that NCBI
      # returns identically twice in a row (e.g. the final chunk, or records counted by
      # esearch but never returned); a truncated transfer gives a different count each time
      if (length(bad) == 0 && (n >= want || (n > 0 && n == prev_n))) {
        writeLines(lines, chunk_file)
        unlink(part)
        ok <- TRUE; got <- n
        break
      }
      reason <- sprintf("got %d of %.0f accessions", n, want)
      if (length(bad)) {
        reason <- sprintf("%s; %d line(s) not accessions: %s", reason, length(bad), substr(bad[1], 1, 120))
        prev_n <- -1
      } else prev_n <- n
    } else {
      reason <- r$msg
    }
    if (try < max_tries) {
      message(sprintf("    %s: attempt %d/%d failed (%s); retrying in %ds",
                      label, try, max_tries, reason, backoff(try)))
      Sys.sleep(backoff(try))
      # every 3rd failure, start a fresh NCBI search session in case the old one expired
      if (try %% 3 == 0) {
        message("    refreshing the NCBI search session...")
        new_srch <- tryCatch(run_esearch(), error = function(e) NULL)
        if (!is.null(new_srch)) {
          if (new_srch$total != total)
            message(sprintf("    WARNING: NCBI now reports %.0f records (was %.0f); the database changed during the download.",
                            new_srch$total, total))
          srch$query_key <- new_srch$query_key
          srch$web_env   <- new_srch$web_env
        }
      }
    }
  }
  unlink(paste0(chunk_file, ".part"))
  if (!ok) stop(sprintf(paste0("%s failed after %d attempts (%s).\n",
                               "  Completed chunks are kept in %s; re-run the same command to resume."),
                        label, max_tries, reason, acc_dir), call. = FALSE)
  got <- min(got, want)   # extra lines are overlap
  if (got < want) {
    message(sprintf("  %s: %d accessions (NCBI consistently returns %.0f fewer than its count; accepted)",
                    label, got, want - got))
  } else message(sprintf("  %s: %d accessions", label, got))
  Sys.sleep(pause)
}

# concatenate chunks into a single accession list (dropping any duplicates)
acc_file <- file.path(acc_dir, paste0(taxon, "_ncbi_acc.txt"))
all_accs <- unique(unlist(lapply(chunk_files, readLines, warn = FALSE)))
writeLines(all_accs, acc_file)
invisible(file.remove(c(chunk_files, state_file)))
message(sprintf("  Accession list: %d unique entries -> %s", length(all_accs), acc_file))
if (length(all_accs) != total)
  message(sprintf(paste0("  Note: NCBI's search counted %.0f records; %d unique accessions were retrieved.\n",
                         "        (NCBI's count includes some records it does not return, e.g. suppressed or\n",
                         "        withdrawn ones; records can also change at NCBI during a long download.)"),
                  total, length(all_accs)))

# =======================================================================
# STEP 3 — build the subset BLAST database
# =======================================================================
message("\n--- Step 3: Building subset BLAST database ---")
# pass the search details and this wrapper's path to the subset script for its readme
wrapper_path <- if (!is.na(script_path)) script_path else ""
Sys.setenv(
  SUBSET_SEARCH_QUERY   = query,
  SUBSET_SEARCH_DATE    = search_date,
  SUBSET_SEARCH_COUNT   = sprintf("%.0f", total),
  SUBSET_CALLER_SCRIPT  = wrapper_path,
  SUBSET_CALLER_COMMAND = paste(c("Rscript", shQuote(wrapper_path), shQuote(args)), collapse = " ")
)
run(
  paste(
    "Rscript", shQuote(subset_script),
    shQuote(db_path),
    shQuote(acc_file),
    shQuote(out_dir),
    shQuote(title)
  ),
  "subset_blastdb_by_accession.R"
)

message("\n=== Done ===")
message("  Accession list: ", acc_file)
message("  readme.txt:     ", file.path(out_dir, "readme.txt"))
message("  Database:       ", file.path(out_dir, gsub("[^A-Za-z0-9._-]", "_", basename(out_dir))), ".*")
