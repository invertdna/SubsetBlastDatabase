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

# ---- ensure edirect is on PATH ----------------------------------------
edirect_bin <- file.path(Sys.getenv("HOME"), "edirect")
Sys.setenv(PATH = paste(Sys.getenv("PATH"), edirect_bin, sep = ":"))

# ---- helper: run a shell command, stop on failure ---------------------
run <- function(cmd, label = cmd) {
  message("  cmd: ", cmd)
  status <- system(cmd)
  if (status != 0) stop(label, " failed with exit code ", status)
  invisible(status)
}

# =======================================================================
# STEP 1 — esearch: get total count and save WebEnv/query_key for paging
# =======================================================================
message("\n--- Step 1: Querying NCBI nuccore ---")
message("  query: ", query)

search_xml <- tempfile(fileext = ".xml")
run(
  paste0("esearch -db nuccore -query ", shQuote(query), " > ", shQuote(search_xml)),
  "esearch"
)
search_date <- format(Sys.time(), "%Y-%m-%d %H:%M %Z")

xml_lines <- readLines(search_xml)
count_line <- grep("<Count>", xml_lines, value = TRUE)
if (length(count_line) == 0) {
  stop("Could not parse record count. esearch output:\n", paste(xml_lines, collapse = "\n"))
}
total <- as.integer(sub(".*<Count>(\\d+)</Count>.*", "\\1", count_line[1]))
message(sprintf("  Total matching records: %d", total))
if (total == 0) stop("No records match the query.")

# =======================================================================
# STEP 2 — efetch accessions in chunks, using stored search history
# =======================================================================
message("\n--- Step 2: Fetching accessions ---")

chunk_size <- 100000
n_chunks   <- ceiling(total / chunk_size)
acc_dir    <- file.path(out_dir, "accessions")
dir.create(acc_dir, showWarnings = FALSE, recursive = TRUE)

chunk_files <- character(n_chunks)
for (i in seq_len(n_chunks)) {
  start <- (i - 1) * chunk_size
  end   <- min(i * chunk_size - 1, total - 1)
  message(sprintf("  Chunk %d/%d  (records %d – %d)...", i, n_chunks, start, end))

  chunk_file <- file.path(acc_dir, sprintf("%s_chunk%03d.txt", taxon, i))
  run(
    paste0(
      "cat ", shQuote(search_xml),
      " | efetch -format acc -start ", start, " -stop ", end,
      " > ", shQuote(chunk_file)
    ),
    paste0("efetch chunk ", i)
  )

  n_saved <- length(readLines(chunk_file))
  message(sprintf("    -> %d accessions", n_saved))
  chunk_files[i] <- chunk_file

  if (i < n_chunks) Sys.sleep(1)   # be polite to NCBI
}

invisible(file.remove(search_xml))

# concatenate chunks into a single accession list
acc_file  <- file.path(acc_dir, paste0(taxon, "_ncbi_acc.txt"))
all_accs  <- unlist(lapply(chunk_files, readLines))
all_accs  <- all_accs[nzchar(trimws(all_accs))]
writeLines(all_accs, acc_file)
invisible(file.remove(chunk_files))
message(sprintf("  Accession list: %d entries -> %s", length(all_accs), acc_file))

# =======================================================================
# STEP 3 — build the subset BLAST database
# =======================================================================
message("\n--- Step 3: Building subset BLAST database ---")
# pass the search details and this wrapper's path to the subset script for its readme
wrapper_path <- if (!is.na(script_path)) script_path else ""
Sys.setenv(
  SUBSET_SEARCH_QUERY   = query,
  SUBSET_SEARCH_DATE    = search_date,
  SUBSET_SEARCH_COUNT   = as.character(total),
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
