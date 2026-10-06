#!/usr/bin/env Rscript
# demo_balaenoptera_hyak.R
# Hyak (UW klone) version of demo_balaenoptera.R. Demonstrates the SubsetBlastDatabase
# pipeline with BLAST+ run from an Apptainer container:
#   1. Build a Balaenoptera subset of core_nt via make_subset_blastdb_hyak.sh
#      (fetches accessions from NCBI, builds the standalone database, and records
#      provenance: readme.txt, code/ copies of the scripts, command.txt, and the
#      accession list in accessions/)
#   2. Show the provenance record (readme.txt) written for the new database
#   3. Search the subset database with testquery.fasta (100 % identity, top 5 hits)
#
# Usage:  Rscript demo_balaenoptera_hyak.R      (no arguments; also works via source())
#
# Does NOT read ~/.bashrc. Settings come from environment variables, which you can
# export in your sbatch script before calling Rscript, or set in the "settings" block below:
#   BLAST_SIF      path to the BLAST+ Apptainer image
#                  (get one with: apptainer pull blast.sif docker://ncbi/blast:latest)
#   NCBI_EMAIL     your email (NCBI asks E-utilities users to identify themselves)
#   NCBI_API_KEY   NCBI API key (recommended)
#
# How BLAST+ is run (same rules as make_subset_blastdb_hyak.sh):
#   - if R itself is running inside a container (APPTAINER_CONTAINER is set), BLAST+
#     must be on that container's PATH and is called directly;
#   - otherwise, if BLAST_SIF is set, BLAST+ is run with `apptainer exec $BLAST_SIF ...`;
#   - otherwise BLAST+ must be on the host PATH (e.g. after `module load`).
#
# On klone, apptainer is only available on compute nodes, so run this in a Slurm job
# (sbatch or salloc), not on a login node. Example sbatch script:
#
#   #!/bin/bash
#   #SBATCH --job-name=demo_balaenoptera
#   #SBATCH --account=<your_account>
#   #SBATCH --partition=<your_partition>
#   #SBATCH --cpus-per-task=8
#   #SBATCH --mem=32G
#   #SBATCH --time=4:00:00
#   export BLAST_SIF=/gscratch/<group>/containers/blast.sif
#   export NCBI_EMAIL="you@uw.edu"
#   export NCBI_API_KEY="your_key"
#   Rscript /path/to/SubsetBlastDatabase/demo_balaenoptera_hyak.R
#
# Edit `blastdb` below to point at core_nt on Hyak.

# ---- settings --------------------------------------------------------------
blastdb   <- "/mmfs1/home/rpkelly/rpkelly/core_nt/core_nt"   # source database on Hyak (path + db name)
# Fill these in here only if you don't export them in your job script
# (values already in the environment take precedence):
settings <- c(
  BLAST_SIF    = "/gscratch/rpkelly/apptainers/kellylab.sif",   # e.g. "/gscratch/<group>/containers/blast.sif"
  NCBI_EMAIL   = "rpkelly@uw.edu",   # e.g. "you@uw.edu"
  NCBI_API_KEY = "213c27165c0e007db4c206d42385970be108"
)
for (k in names(settings)) {
  if (nzchar(settings[[k]]) && !nzchar(Sys.getenv(k))) do.call(Sys.setenv, setNames(list(settings[[k]]), k))
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
taxon        <- "Balaenoptera"
query        <- "(Balaenoptera[Organism]) AND (mitochondrion OR mitochondrial) AND 100:20000[Sequence Length]"
title        <- "Balaenoptera mitochondrial DNA"
query_fasta  <- file.path(script_dir, "testquery.fasta")
out_dir      <- file.path(script_dir, taxon)
blast_out    <- file.path(out_dir, "testquery_blast_results.txt")
# database files are named after the output folder, with BLAST-unsafe characters -> "_"
db_name      <- gsub("[^A-Za-z0-9._-]", "_", basename(out_dir))
make_script  <- file.path(script_dir, "make_subset_blastdb_hyak.sh")

if (!dir.exists(dirname(blastdb)))
  stop("Source database folder not found: ", dirname(blastdb),
       "\n  (Edit `blastdb` at the top of this script.)")
if (!file.exists(query_fasta)) stop("Query file not found: ", query_fasta)
if (!file.exists(make_script)) stop("Cannot find make_subset_blastdb_hyak.sh in: ", script_dir)

# taxdb files live alongside core_nt; add that directory to BLASTDB
# (passed through to the container, since apptainer keeps the host environment)
Sys.setenv(BLASTDB = paste(normalizePath(dirname(blastdb)), Sys.getenv("BLASTDB"), sep = ":"))

# ---- how to run BLAST+ (for step 3; the shell script works this out itself) ----
in_container <- nzchar(Sys.getenv("APPTAINER_CONTAINER")) || nzchar(Sys.getenv("SINGULARITY_CONTAINER"))
blast_sif    <- Sys.getenv("BLAST_SIF")
if (!in_container && nzchar(blast_sif)) {
  if (!file.exists(blast_sif)) stop("BLAST_SIF is set but the image does not exist: ", blast_sif)
  apptainer <- Sys.getenv("APPTAINER_BIN")
  if (!nzchar(apptainer)) apptainer <- Sys.which("apptainer")
  if (!nzchar(apptainer)) apptainer <- Sys.which("singularity")
  if (!nzchar(apptainer))
    stop("apptainer not found. On klone it is only available on compute nodes;\n",
         "  run this inside a Slurm job (sbatch or salloc), not on a login node.")
  blast_mode <- paste(apptainer, "exec", blast_sif)
} else if (nzchar(Sys.which("blastn"))) {
  blast_mode <- if (in_container) "inside container" else "host PATH"
} else {
  stop("BLAST+ not found. Set BLAST_SIF to a BLAST+ Apptainer image ",
       "(apptainer pull blast.sif docker://ncbi/blast:latest),\n",
       "  or run R inside a container that provides BLAST+.")
}
message("BLAST+ will be run via: ", blast_mode)

# ---- step 1: build Balaenoptera subset database -------------------------
message("=== Step 1: Building ", taxon, " subset database ===")
status <- system(paste(
  "bash", shQuote(make_script),
  shQuote(blastdb),
  shQuote(taxon),
  shQuote(query),
  shQuote(title),
  shQuote(out_dir)
))
if (status != 0) stop("make_subset_blastdb_hyak.sh failed with exit code ", status)

# ---- step 2: show the provenance record ----------------------------------
message("\n=== Step 2: Provenance record (", file.path(taxon, "readme.txt"), ") ===")
writeLines(readLines(file.path(out_dir, "readme.txt")))
message("\nCode and exact commands saved in: ", file.path(out_dir, "code"))
message(paste0("  ", list.files(file.path(out_dir, "code")), collapse = "\n"))

# ---- step 3: blastn search (100 % identity, max 5 hits) -----------------
message("\n=== Step 3: Running blastn on subset database ===")
# BLAST+ splits -db names on spaces, so run blastn from inside the database folder
# and use the bare database name.
blastn_args <- paste(
  "-db",          shQuote(db_name),
  "-query",       shQuote(query_fasta),
  "-perc_identity 100",
  "-max_target_seqs 5",
  "-outfmt",      shQuote("6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"),
  "-out",         shQuote(blast_out)
)
if (!in_container && nzchar(blast_sif)) {
  # make the database folder, the query file's folder and the core_nt folder visible
  # inside the container
  binds <- unique(normalizePath(c(out_dir, dirname(query_fasta), dirname(blastdb))))
  blast_cmd <- paste(
    "cd", shQuote(out_dir), "&&",
    shQuote(apptainer), "exec", "--bind", shQuote(paste(binds, collapse = ",")),
    shQuote(blast_sif), "blastn", blastn_args
  )
} else {
  blast_cmd <- paste("cd", shQuote(out_dir), "&&", "blastn", blastn_args)
}
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
