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

# ---- 1. extract clean fasta: accession.version as sole header ----
#   %a = accession.version, %t = title, %s = sequence
#   -entry_batch accepts accession or accession.version strings directly.
cmd_fasta <- paste0(
  "blastdbcmd",
  " -db ", shQuote(db_path),
  " -entry_batch ", shQuote(acc_file),
  ' -outfmt ">%a %t\n%s"',
  " -out ", shQuote(fasta_file)
)
message("Extracting FASTA sequences...")
message("  cmd: ", cmd_fasta)
status <- system(cmd_fasta)
# exit code 1 is normal when some accessions are absent from the local db (Skipped warnings);
# treat as fatal only if the output file is missing or empty.
if (status != 0) {
  if (!file.exists(fasta_file) || file.size(fasta_file) == 0)
    stop("blastdbcmd (fasta) produced no output (exit code ", status, ")")
  message("  Note: blastdbcmd skipped some accessions not present in local db (exit code ", status, ")")
}

# ---- 2. extract accession and taxid for same entries ----
cmd_tax <- paste0(
  "blastdbcmd",
  " -db ", shQuote(db_path),
  " -entry_batch ", shQuote(acc_file),
  ' -outfmt "%a\t%T"',
  " -out ", shQuote(taxid_file)
)
message("Extracting taxon IDs...")
message("  cmd: ", cmd_tax)
status <- system(cmd_tax)
if (status != 0) {
  if (!file.exists(taxid_file) || file.size(taxid_file) == 0)
    stop("blastdbcmd (taxid) produced no output (exit code ", status, ")")
  message("  Note: blastdbcmd skipped some accessions not present in local db (exit code ", status, ")")
}

# ---- 3. build accession-to-taxid map, deduplicated ----
tax <- read.delim(taxid_file, header = FALSE, colClasses = "character",
                  col.names = c("acc", "taxid"))
tax <- tax[!is.na(tax$taxid) & tax$taxid != "0", ]
tax <- tax[!duplicated(tax$acc), ]
write.table(tax, taxid_map, sep = "\t", row.names = FALSE,
            col.names = FALSE, quote = FALSE)
message(sprintf("Taxid map written: %d entries", nrow(tax)))

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

# ---- 6. clean up intermediate files ----
invisible(file.remove(fasta_file))
message("Removed intermediate FASTA: ", fasta_file)

message("\nDone. New database: ", new_db_name)
message("Files created:")
message("  Taxid list: ", taxid_file)
message("  Taxid map:  ", taxid_map)
message("  BLAST db:   ", new_db_name, ".*")
