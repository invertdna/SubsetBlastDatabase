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

# ---- parse args ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript subset_blastdb_by_accession.R <blastdb> <acc_list> [output_dir] [title]")
}

db_path  <- args[1]
acc_file <- args[2]
out_dir  <- if (length(args) >= 3) args[3] else "subset_db"
db_title <- if (length(args) >= 4) args[4] else basename(out_dir)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# use the output directory name as the file prefix (e.g. "Sebastes" -> "Sebastes.fasta")
prefix      <- basename(out_dir)
fasta_file  <- file.path(out_dir, paste0(prefix, ".fasta"))
taxid_file  <- file.path(out_dir, paste0(prefix, "_taxids.txt"))
taxid_map   <- file.path(out_dir, paste0(prefix, "_taxid_map.txt"))
new_db_name <- file.path(out_dir, prefix)

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

# ---- detect available cores ----
n_cores <- max(1L, parallel::detectCores(logical = TRUE))
message(sprintf("Parallel workers: %d", n_cores))

# ---- split accession list into N chunks, write chunk files ----
n_chunks   <- min(n_cores, length(acc_lines))
chunk_size <- ceiling(length(acc_lines) / n_chunks)
acc_chunk_list <- split(acc_lines, (seq_along(acc_lines) - 1L) %/% chunk_size)
n_chunks   <- length(acc_chunk_list)  # recompute: may be < n_cores for tiny lists

acc_chunk_files <- lapply(seq_len(n_chunks), function(i) {
  f <- file.path(out_dir, sprintf(".%s_acc_chunk_%02d.txt", prefix, i))
  writeLines(acc_chunk_list[[i]], f)
  f
})
message(sprintf("Split %d accessions into %d chunk(s) of up to %d",
                length(acc_lines), n_chunks, chunk_size))

# ---- 1. extract FASTA in parallel ----
message(sprintf("Extracting FASTA sequences (%d parallel job(s))...", n_chunks))
fasta_chunk_files <- lapply(seq_len(n_chunks), function(i)
  file.path(out_dir, sprintf(".%s_fasta_chunk_%02d.fasta", prefix, i)))
# pre-create empty placeholders so cat succeeds even if blastdbcmd finds nothing
invisible(lapply(fasta_chunk_files, file.create))
invisible(parallel::mclapply(seq_len(n_chunks), function(i) {
  cmd <- paste0(
    "blastdbcmd",
    " -db ", shQuote(db_path),
    " -entry_batch ", shQuote(acc_chunk_files[[i]]),
    ' -outfmt ">%a %t\n%s"',
    " -out ", shQuote(fasta_chunk_files[[i]])
  )
  system(cmd)
}, mc.cores = n_chunks))
system(paste("cat", paste(shQuote(unlist(fasta_chunk_files)), collapse = " "),
             ">", shQuote(fasta_file)))
invisible(file.remove(unlist(fasta_chunk_files)))
if (!file.exists(fasta_file) || file.size(fasta_file) == 0)
  stop("blastdbcmd (fasta) produced no output")
message("  FASTA chunks merged")

# ---- 2. extract accession and taxid in parallel ----
message(sprintf("Extracting taxon IDs (%d parallel job(s))...", n_chunks))
taxid_chunk_files <- lapply(seq_len(n_chunks), function(i)
  file.path(out_dir, sprintf(".%s_taxid_chunk_%02d.txt", prefix, i)))
invisible(lapply(taxid_chunk_files, file.create))
invisible(parallel::mclapply(seq_len(n_chunks), function(i) {
  cmd <- paste0(
    "blastdbcmd",
    " -db ", shQuote(db_path),
    " -entry_batch ", shQuote(acc_chunk_files[[i]]),
    ' -outfmt "%a\t%T"',
    " -out ", shQuote(taxid_chunk_files[[i]])
  )
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
cmd_make <- paste0(
  "makeblastdb",
  " -in ", shQuote(fasta_file),
  " -dbtype nucl",
  " -parse_seqids",
  " -taxid_map ", shQuote(taxid_map),
  " -out ", shQuote(new_db_name),
  " -title ", shQuote(db_title)
)
message("Building new BLAST database (title=", db_title, ")...")
message("  cmd: ", cmd_make)
status <- system(cmd_make)
if (status != 0) stop("makeblastdb failed with exit code ", status)

# ---- 6. write readme ----
script_name <- tryCatch(
  basename(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])),
  error = function(e) "subset_blastdb_by_accession.R"
)
readme_file <- file.path(out_dir, "readme.txt")
writeLines(c(
  "SubsetBlastDatabase",
  "===================",
  sprintf("Script:            %s", script_name),
  sprintf("Source database:   %s", db_path),
  sprintf("Created by:        %s", Sys.getenv("USER")),
  sprintf("Created:           %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("Unique accessions: %d", n_seqs),
  sprintf("Unique taxon IDs:  %d", n_taxids)
), readme_file)
message("readme.txt written: ", readme_file)

# ---- 7. clean up intermediate files ----
invisible(file.remove(fasta_file))
message("Removed intermediate FASTA: ", fasta_file)

message("\nDone. New database: ", new_db_name)
message("Files created:")
message("  readme.txt: ", readme_file)
message("  Taxid list: ", taxid_file)
message("  Taxid map:  ", taxid_map)
message("  BLAST db:   ", new_db_name, ".*")
