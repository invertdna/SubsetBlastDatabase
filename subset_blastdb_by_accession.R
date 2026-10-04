#!/usr/bin/env Rscript
# subset_blastdb_by_accession.R
# Create a standalone nucleotide BLAST database from a subset of NCBI accession numbers.
#
# Prefer this over the GI-based version: GI numbers are deprecated by NCBI
# and no longer assigned to new records. Accession.version strings are stable.
#
# Usage:
#   Rscript subset_blastdb_by_accession.R <path_to_blastdb> <path_to_acc_list> [output_dir] [title]
#
# Arguments:
#   path_to_blastdb : path (with db name prefix) to an existing local BLAST db
#                     e.g. /Volumes/Clupea/core_nt/core_nt
#   path_to_acc_list: text file with one accession (or accession.version) per line
#                     e.g. NC_001234.1  or  NC_001234
#   output_dir      : (optional) directory for output files; default "subset_db"
#   title           : (optional) human-readable title for the database; default = output_dir name
#
# The accession list is split into N chunks and blastdbcmd is run in parallel,
# where N = number of logical CPU cores. Uses base R parallel package (no extra deps).
#
# Provenance: a copy of this script, the exact command used, and the input
# accession list are saved to <output_dir>/code/, and a short readme.txt
# describing the database contents and dates is written to <output_dir>/.
# When called from make_subset_blastdb.R, the NCBI search query, date and
# record count (passed via SUBSET_SEARCH_* env vars) are recorded as well,
# and the wrapper script is copied into code/ too.

# ---- load the user's shell PATH from ~/.bashrc ----
# RStudio (including its Terminal and R's system() calls) does not inherit the PATH
# set up for the macOS Terminal app, so BLAST+ and edirect may not be found. Source
# ~/.bashrc in bash and use the resulting PATH for every command this script runs.
if (file.exists(path.expand("~/.bashrc"))) {
  bashrc_path <- suppressWarnings(system(
    "bash -c 'source ~/.bashrc >/dev/null 2>&1; printf %s \"$PATH\"'", intern = TRUE))
  if (length(bashrc_path) == 1 && nzchar(bashrc_path)) Sys.setenv(PATH = bashrc_path)
}

# ---- parse args ----
args <- commandArgs(trailingOnly = TRUE)

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

# full path to this script (for the provenance copy saved with the database)
script_path <- this_script_path()
run_dir     <- getwd()
if (length(args) < 2) {
  stop("Usage: Rscript subset_blastdb_by_accession.R <blastdb> <acc_list> [output_dir] [title]")
}

db_path  <- args[1]
acc_file <- args[2]
out_dir  <- if (length(args) >= 3) args[3] else "subset_db"
db_title <- if (length(args) >= 4) args[4] else basename(out_dir)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_abs <- normalizePath(out_dir)   # absolute path; used for all intermediate files

# ---- BLAST+ and paths with spaces ----
# BLAST+ splits database names on spaces (-db "a b" means two databases), so a
# database path under e.g. "Kelly_Lab - Documents/2. KellyLab" breaks blastdbcmd
# and makeblastdb. Workaround: run each BLAST command from inside the database's
# own folder and refer to the database by its bare (space-free) name.
in_dir  <- function(dir, cmd) paste("cd", shQuote(dir), "&&", cmd)
db_dir  <- normalizePath(dirname(db_path), mustWork = TRUE)
db_name <- basename(db_path)
if (grepl("[[:space:]]", db_name))
  stop("The source database name itself contains spaces, which BLAST+ cannot handle: ", db_name)

# use the output directory name as the file prefix (e.g. "Sebastes" -> "Sebastes.fasta");
# characters BLAST+ can't handle in a db name (spaces etc.) become "_"
prefix      <- gsub("[^A-Za-z0-9._-]", "_", basename(out_dir))
fasta_file  <- file.path(out_abs, paste0(prefix, ".fasta"))
taxid_file  <- file.path(out_abs, paste0(prefix, "_taxids.txt"))
taxid_map   <- file.path(out_abs, paste0(prefix, "_taxid_map.txt"))
new_db_name <- file.path(out_abs, prefix)

# ---- 0. validate input: warn if entries look like raw GI numbers ----
#   GI numbers are pure integers. Accessions contain letters (e.g. NC_001234.1).
acc_lines <- readLines(acc_file)
acc_lines <- acc_lines[nzchar(trimws(acc_lines))]  # drop blank lines
gi_like   <- grepl("^[0-9]+$", trimws(acc_lines))
if (any(gi_like)) {
  warning(sprintf(
    "%d line(s) look like GI numbers (pure integers). ",
    sum(gi_like),
    "Use subset_blastdb_by_GI.R for GI-based lookups, or verify your accession list."
  ))
}
message(sprintf("Input: %d accession(s) from %s", length(acc_lines), acc_file))

# ---- check the source database can be opened before doing any work ----
db_info <- suppressWarnings(
  system(in_dir(db_dir, paste("blastdbcmd -db", shQuote(db_name), "-info 2>&1")), intern = TRUE))
if (!is.null(attr(db_info, "status")) && attr(db_info, "status") != 0) {
  # remove the output dir only if it is empty (i.e. we just created it)
  if (length(list.files(out_abs, all.files = TRUE, no.. = TRUE)) == 0) unlink(out_abs, recursive = TRUE)
  stop("Cannot open the source BLAST database: ", db_path, "\n  blastdbcmd said: ",
       paste(head(db_info, 3), collapse = "\n    "),
       "\n  Common causes: the drive holding the database is reached over a network share\n  (SMB/NFS -- BLAST's LMDB index needs a locally attached disk), the download of the\n  database is incomplete, or the path is wrong.", call. = FALSE)
}
message("Source database opened OK: ", db_path)

# ---- detect available cores ----
n_cores <- max(1L, parallel::detectCores(logical = TRUE))
message(sprintf("Parallel workers: %d", n_cores))

# ---- split accession list into N chunks, write chunk files ----
n_chunks   <- min(n_cores, length(acc_lines))
chunk_size <- ceiling(length(acc_lines) / n_chunks)
acc_chunk_list <- split(acc_lines, (seq_along(acc_lines) - 1L) %/% chunk_size)
n_chunks   <- length(acc_chunk_list)  # recompute: may be < n_cores for tiny lists

acc_chunk_files <- lapply(seq_len(n_chunks), function(i) {
  f <- file.path(out_abs, sprintf(".%s_acc_chunk_%02d.txt", prefix, i))
  writeLines(acc_chunk_list[[i]], f)
  f
})
message(sprintf("Split %d accessions into %d chunk(s) of up to %d",
                length(acc_lines), n_chunks, chunk_size))

# ---- 1. extract FASTA in parallel ----
message(sprintf("Extracting FASTA sequences (%d parallel job(s))...", n_chunks))
fasta_chunk_files <- lapply(seq_len(n_chunks), function(i)
  file.path(out_abs, sprintf(".%s_fasta_chunk_%02d.fasta", prefix, i)))
# pre-create empty placeholders so cat succeeds even if blastdbcmd finds nothing
invisible(lapply(fasta_chunk_files, file.create))
invisible(parallel::mclapply(seq_len(n_chunks), function(i) {
  cmd <- in_dir(db_dir, paste0(
    "blastdbcmd",
    " -db ", shQuote(db_name),
    " -entry_batch ", shQuote(acc_chunk_files[[i]]),
    ' -outfmt ">%a %t\n%s"',
    " -out ", shQuote(fasta_chunk_files[[i]])
  ))
  system(cmd)
}, mc.cores = n_chunks))
system(paste("cat", paste(shQuote(unlist(fasta_chunk_files)), collapse = " "),
             ">", shQuote(fasta_file)))
invisible(file.remove(unlist(fasta_chunk_files)))
if (!file.exists(fasta_file) || file.size(fasta_file) == 0) {
  invisible(file.remove(c(fasta_file, unlist(acc_chunk_files))))
  stop("blastdbcmd (fasta) produced no output: none of the accessions were found in ", db_path)
}
message("  FASTA chunks merged")

# ---- 2. extract accession and taxid in parallel ----
message(sprintf("Extracting taxon IDs (%d parallel job(s))...", n_chunks))
taxid_chunk_files <- lapply(seq_len(n_chunks), function(i)
  file.path(out_abs, sprintf(".%s_taxid_chunk_%02d.txt", prefix, i)))
invisible(lapply(taxid_chunk_files, file.create))
invisible(parallel::mclapply(seq_len(n_chunks), function(i) {
  cmd <- in_dir(db_dir, paste0(
    "blastdbcmd",
    " -db ", shQuote(db_name),
    " -entry_batch ", shQuote(acc_chunk_files[[i]]),
    ' -outfmt "%a\t%T"',
    " -out ", shQuote(taxid_chunk_files[[i]])
  ))
  system(cmd)
}, mc.cores = n_chunks))
system(paste("cat", paste(shQuote(unlist(taxid_chunk_files)), collapse = " "),
             ">", shQuote(taxid_file)))
invisible(file.remove(c(unlist(taxid_chunk_files), unlist(acc_chunk_files))))
if (!file.exists(taxid_file) || file.size(taxid_file) == 0)
  stop("blastdbcmd (taxid) produced no output")
message("  Taxid chunks merged")

# ---- 3. build accession-to-taxid map, deduplicated ----
tax <- read.delim(taxid_file, header = FALSE, colClasses = "character",
                  col.names = c("acc", "taxid"))
tax <- tax[!is.na(tax$taxid) & tax$taxid != "0", ]
tax <- tax[!duplicated(tax$acc), ]
write.table(tax, taxid_map, sep = "\t", row.names = FALSE,
            col.names = FALSE, quote = FALSE)
message(sprintf("Taxid map written: %d entries", nrow(tax)))
n_taxids <- length(unique(tax$taxid))

# ---- 4. deduplicate fasta ----
#   Multiple input accessions may expand to the same sequence; makeblastdb
#   -parse_seqids will error on duplicate seq IDs, so keep only first occurrence.
lines  <- readLines(fasta_file)
is_hdr <- grepl("^>", lines)
ids    <- sub("^>(\\S+).*", "\\1", lines[is_hdr])
grp    <- cumsum(is_hdr)
keep   <- !duplicated(ids)
keep_idx  <- which(is_hdr)[keep]
keep_mask <- grp %in% grp[keep_idx]
writeLines(lines[keep_mask], fasta_file)
n_seqs  <- sum(keep)
n_dedup <- sum(!keep)
if (n_dedup > 0) message(sprintf("Removed %d duplicate accession(s) from FASTA", n_dedup))

# ---- 5. build new standalone blast database ----
# run inside out_abs with bare file names (see "BLAST+ and paths with spaces" above)
cmd_make <- in_dir(out_abs, paste0(
  "makeblastdb",
  " -in ", shQuote(basename(fasta_file)),
  " -dbtype nucl",
  " -parse_seqids",
  " -taxid_map ", shQuote(basename(taxid_map)),
  " -out ", shQuote(prefix),
  " -title ", shQuote(db_title)
))
message("Building new BLAST database (title=", db_title, ")...")
message("  cmd: ", cmd_make)
status <- system(cmd_make)
if (status != 0) stop("makeblastdb failed with exit code ", status)

# ---- 6. save code + write readme ----
run_user <- Sys.info()[["user"]]
# copy of the generating code, the exact command, and the input accession list
code_dir <- file.path(out_dir, "code")
dir.create(code_dir, showWarnings = FALSE, recursive = TRUE)
script_name <- if (!is.na(script_path)) basename(script_path) else "subset_blastdb_by_accession.R"
if (!is.na(script_path) && file.exists(script_path)) {
  invisible(file.copy(script_path, file.path(code_dir, script_name), overwrite = TRUE))
} else {
  warning("Could not locate this script's path (run via Rscript to save a copy of the code).")
}
code_line <- paste0("code/", script_name)
caller_script <- Sys.getenv("SUBSET_CALLER_SCRIPT")
if (nzchar(caller_script) && file.exists(caller_script)) {
  invisible(file.copy(caller_script, file.path(code_dir, basename(caller_script)), overwrite = TRUE))
  code_line <- sprintf("code/%s (calls %s)", basename(caller_script), code_line)
}

# keep the accession list: reference it if already inside out_dir, else copy it
acc_abs <- normalizePath(acc_file)
out_abs <- normalizePath(out_dir)
if (startsWith(acc_abs, paste0(out_abs, "/"))) {
  acc_rel <- substring(acc_abs, nchar(out_abs) + 2)
} else {
  invisible(file.copy(acc_file, file.path(code_dir, "input_accessions.txt"), overwrite = TRUE))
  acc_rel <- sprintf("code/input_accessions.txt (copy of %s)", basename(acc_file))
}

blast_version <- tryCatch(
  system("blastdbcmd -version", intern = TRUE)[1],
  error = function(e) "unknown", warning = function(w) "unknown"
)
writeLines(c(
  sprintf("# Run from: %s", run_dir),
  sprintf("# Run on:   %s by %s@%s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
          run_user, Sys.info()[["nodename"]]),
  sprintf("# BLAST+:   %s", blast_version),
  sprintf("# R:        %s", R.version.string),
  if (nzchar(Sys.getenv("SUBSET_CALLER_COMMAND")))
    c("# Top-level command:", Sys.getenv("SUBSET_CALLER_COMMAND"), "# which ran:"),
  paste(c("Rscript", shQuote(if (is.na(script_path)) script_name else script_path),
          shQuote(args)), collapse = " ")
), file.path(code_dir, "command.txt"))

# date of the source database snapshot (from blastdbcmd -info)
# (db_info was read by the upfront check at the start of the script)
date_line   <- grep("Date:", db_info, value = TRUE)[1]
src_db_date <- if (is.na(date_line)) "unknown" else
  trimws(sub(".*Date:\\s*(.*[AP]M).*", "\\1", date_line))

readme_file <- file.path(out_dir, "readme.txt")
writeLines(c(
  db_title,
  sprintf("Standalone nucleotide BLAST database (%s.*), built %s by %s.",
          prefix, format(Sys.Date(), "%Y-%m-%d"), run_user),
  sprintf("Contents: %d unique accessions, %d unique NCBI taxon IDs,", n_seqs, n_taxids),
  sprintf("  from the accession list %s.", acc_rel),
  if (nzchar(Sys.getenv("SUBSET_SEARCH_QUERY"))) c(
    sprintf("Search: NCBI nuccore, run %s, %s records matched:",
            Sys.getenv("SUBSET_SEARCH_DATE", "unknown date"),
            Sys.getenv("SUBSET_SEARCH_COUNT", "?")),
    paste0("  ", Sys.getenv("SUBSET_SEARCH_QUERY"))),
  sprintf("Source: %s (source db dated %s).", db_path, src_db_date),
  sprintf("Code: %s; exact command(s) in code/command.txt.", code_line)
), readme_file)
message("readme.txt written: ", readme_file)

# ---- 7. clean up intermediate files ----
invisible(file.remove(fasta_file))
message("Removed intermediate FASTA: ", fasta_file)

message("\nDone. New database: ", new_db_name)
message("Files created:")
message("  readme.txt: ", readme_file)
message("  Code:       ", code_dir, "/")
message("  Taxid list: ", taxid_file)
message("  Taxid map:  ", taxid_map)
message("  BLAST db:   ", new_db_name, ".*")
