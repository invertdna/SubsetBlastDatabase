#!/usr/bin/env Rscript
# demo_balaenoptera.R
# Demonstrates the SubsetBlastDatabase pipeline:
#   1. Build a Balaenoptera subset of core_nt via make_subset_blastdb.R
#   2. Search the subset database with testquery.fasta (100 % identity, top 5 hits)

script_dir <- dirname(normalizePath(
  sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])
))

# ---- paths ---------------------------------------------------------------
blastdb      <- "/Volumes/Clupea/core_nt/core_nt"
query_fasta  <- file.path(script_dir, "testquery.fasta")
out_dir      <- file.path(script_dir, "Balaenoptera")
blast_out    <- file.path(out_dir, "testquery_blast_results.txt")
db_prefix    <- file.path(out_dir, "Balaenoptera")

# taxdb files live alongside core_nt; add that directory to BLASTDB
Sys.setenv(BLASTDB = paste("/Volumes/Clupea/core_nt", Sys.getenv("BLASTDB"), sep = ":"))

# ---- step 1: build Balaenoptera subset database -------------------------
message("=== Step 1: Building Balaenoptera subset database ===")
status <- system(paste(
  "Rscript", shQuote(file.path(script_dir, "make_subset_blastdb.R")),
  shQuote(blastdb),
  "Balaenoptera",
  shQuote("(Balaenoptera[Organism]) AND (mitochondrion OR mitochondrial) AND 100:20000[Sequence Length]"),
  shQuote("Balaenoptera mitochondrial DNA"),
  shQuote(out_dir)
))
if (status != 0) stop("make_subset_blastdb.R failed with exit code ", status)

# ---- step 2: blastn search (100 % identity, max 5 hits) -----------------
message("\n=== Step 2: Running blastn on subset database ===")
blast_cmd <- paste(
  "blastn",
  "-db",          shQuote(db_prefix),
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
results <- read.table(blast_out, sep = "\t", col.names = cols)
if (nrow(results) == 0) {
  message("No hits at 100 % identity.")
} else {
  print(results)
  message(sprintf("\n%d hit(s) written to: %s", nrow(results), blast_out))
}
