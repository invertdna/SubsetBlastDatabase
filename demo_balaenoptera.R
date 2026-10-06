#!/usr/bin/env Rscript
# demo_balaenoptera.R
# Demonstrates the SubsetBlastDatabase pipeline:
#   1. Build a Balaenoptera subset of core_nt via make_subset_blastdb.R
#      (fetches accessions from NCBI, builds the standalone database, and records
#      provenance: readme.txt, code/ copies of the scripts, command.txt, and the
#      accession list in accessions/)
#   2. Show the provenance record (readme.txt) written for the new database
#   3. Search the subset database with testquery.fasta (100 % identity, top 5 hits)
#
# Usage:  Rscript demo_balaenoptera.R      (no arguments; also works via source())
#
# Paths are resolved relative to this script's folder, and are safe for OneDrive-style
# paths with spaces (e.g. "Kelly_Lab - Documents/2. KellyLab"). Edit `blastdb` below
# to point at your local copy of core_nt.

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

script_path <- this_script_path()
if (is.na(script_path)) stop("Cannot determine this script's location; run with Rscript or source().")
script_dir <- dirname(script_path)

# ---- paths ---------------------------------------------------------------
blastdb      <- "/Volumes/Clupea/core_nt/core_nt"   # local source database (path + db name)
taxon        <- "Balaenoptera"
query        <- "(Balaenoptera[Organism]) AND (mitochondrion OR mitochondrial) AND 100:20000[Sequence Length]"
title        <- "Balaenoptera mitochondrial DNA"
query_fasta  <- file.path(script_dir, "testquery.fasta")
out_dir      <- file.path(script_dir, taxon)
blast_out    <- file.path(out_dir, "testquery_blast_results.txt")
# database files are named after the output folder, with BLAST-unsafe characters -> "_"
db_name      <- gsub("[^A-Za-z0-9._-]", "_", basename(out_dir))

if (!dir.exists(dirname(blastdb)))
  stop("Source database folder not found: ", dirname(blastdb),
       "\n  (Is the drive mounted? Edit `blastdb` at the top of this script.)")
if (!file.exists(query_fasta)) stop("Query file not found: ", query_fasta)

# taxdb files live alongside core_nt; add that directory to BLASTDB
Sys.setenv(BLASTDB = paste(dirname(blastdb), Sys.getenv("BLASTDB"), sep = ":"))

# ---- step 1: build Balaenoptera subset database -------------------------
message("=== Step 1: Building ", taxon, " subset database ===")
status <- system(paste(
  "Rscript", shQuote(file.path(script_dir, "make_subset_blastdb.R")),
  shQuote(blastdb),
  shQuote(taxon),
  shQuote(query),
  shQuote(title),
  shQuote(out_dir)
))
if (status != 0) stop("make_subset_blastdb.R failed with exit code ", status)

# ---- step 2: show the provenance record ----------------------------------
message("\n=== Step 2: Provenance record (", file.path(taxon, "readme.txt"), ") ===")
writeLines(readLines(file.path(out_dir, "readme.txt")))
message("\nCode and exact commands saved in: ", file.path(out_dir, "code"))
message(paste0("  ", list.files(file.path(out_dir, "code")), collapse = "\n"))

# ---- step 3: blastn search (100 % identity, max 5 hits) -----------------
message("\n=== Step 3: Running blastn on subset database ===")
# BLAST+ splits -db names on spaces, so run blastn from inside the database folder
# and use the bare database name (the OneDrive path contains spaces)
blast_cmd <- paste(
  "cd", shQuote(out_dir), "&&",
  "blastn",
  "-db",          shQuote(db_name),
  "-query",       shQuote(query_fasta),
  "-perc_identity 100",
  "-max_target_seqs 5",
  "-outfmt",      shQuote("6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"),
  "-out",         shQuote(blast_out)
)
message("  cmd: ", blast_cmd)
status <- system(blast_cmd)
if (status != 0) stop("blastn failed with exit code ", status)

# ---- display results -----------------------------------------------------
message("\n=== BLAST results (100 % identity, max 5 hits per query) ===")
cols <- c("qseqid", "sseqid", "pident", "length", "mismatch",
          "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
if (!file.exists(blast_out) || file.size(blast_out) == 0) {
  message("No hits at 100 % identity.")
} else {
  results <- read.table(blast_out, sep = "\t", col.names = cols)
  print(results)
  message(sprintf("\n%d hit(s) written to: %s", nrow(results), blast_out))
}
