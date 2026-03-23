# SubsetBlastDatabase

Create portable, standalone local BLAST databases from taxonomic subsets of a larger
database (e.g. NCBI `core_nt`).

---

## Scripts

Each script has an R version (`.R`) and a bash equivalent (`.sh`) with identical interfaces.
Use the `.sh` versions on systems where R is not available.

### `make_subset_blastdb.R` / `make_subset_blastdb.sh` — full pipeline wrapper

Fetches matching accession numbers from NCBI and builds the subset database in one step.

**Usage**

```bash
Rscript make_subset_blastdb.R <blastdb> <taxon> <query> <title> [output_dir]
```

| Argument | Description |
|---|---|
| `blastdb` | Path (with db prefix) to an existing local BLAST database |
| `taxon` | Short name used for output file naming — no spaces (e.g. `Sebastes`) |
| `query` | NCBI `esearch` query string |
| `title` | Human-readable title embedded in the output database metadata |
| `output_dir` | (optional) Output directory; defaults to the `taxon` value |

**Example**

```bash
# R version
Rscript make_subset_blastdb.R \
  /Volumes/Clupea/core_nt/core_nt \
  Sebastes \
  "(Sebastes[Organism]) AND (mitochondrion OR mitochondrial) AND 100:20000[Sequence Length]" \
  "Sebastes mitochondrial DNA"

# bash version (no R required)
bash make_subset_blastdb.sh \
  /Volumes/Clupea/core_nt/core_nt \
  Sebastes \
  "(Sebastes[Organism]) AND (mitochondrion OR mitochondrial) AND 100:20000[Sequence Length]" \
  "Sebastes mitochondrial DNA"
```

**What it does**

1. Runs `esearch` once against NCBI nuccore and saves the search history (WebEnv/query_key).
2. Pages through results with `efetch -format acc` in chunks of 100,000, with a 1-second
   pause between chunks to respect NCBI rate limits.
3. Saves the combined accession list to `<output_dir>/accessions/<taxon>_ncbi_acc.txt`.
4. Calls `subset_blastdb_by_accession.R` (R version) or `subset_blastdb_by_accession.sh`
   (bash version) to extract sequences from the local database and build the standalone
   BLAST database.

**Requirements:** see Dependencies below.

---

### `subset_blastdb_by_accession.R` / `subset_blastdb_by_accession.sh` — database builder

Extracts sequences for a list of accession numbers from an existing local BLAST database
and builds a new, self-contained database. Can be used independently of the wrapper.

**Usage**

```bash
Rscript subset_blastdb_by_accession.R <blastdb> <acc_list> [output_dir] [title]
```

| Argument | Description |
|---|---|
| `blastdb` | Path (with db prefix) to an existing local BLAST database |
| `acc_list` | Text file with one accession or accession.version per line |
| `output_dir` | (optional) Output directory; defaults to `subset_db` |
| `title` | (optional) Database title; defaults to the output directory name |

**What it does**

1. Validates the accession list (warns if any lines look like GI numbers).
2. Extracts FASTA sequences from the local database with `blastdbcmd -entry_batch`.
3. Extracts accession-to-taxid mappings with a second `blastdbcmd` call.
4. Deduplicates the FASTA (multi-volume databases return the same sequence once per
   volume; `makeblastdb` requires unique sequence IDs).
5. Builds the new database with `makeblastdb -parse_seqids -taxid_map`.
6. Removes the intermediate FASTA file.

Accessions not present in the local database are silently skipped by `blastdbcmd`
(this is expected when the local database snapshot is older than the NCBI query).

**Output files** (all prefixed with the output directory name):

```
<output_dir>/
  <taxon>_taxids.txt       raw accession-taxid pairs from blastdbcmd
  <taxon>_taxid_map.txt    deduplicated map used by makeblastdb
  <taxon>.{nhr,nin,...}    BLAST database files
```

---

## Dependencies

| Tool | Used by | Notes |
|---|---|---|
| Tool | Used by | Notes |
|---|---|---|
| **R** (≥ 4.0) | `.R` scripts only | base R only; no packages required — not needed if using `.sh` versions |
| **bash** (≥ 4.0), **awk**, **grep** | `.sh` scripts | standard on macOS/Linux; macOS ships bash 3.2 — install bash 4+ via Homebrew if needed |
| **NCBI BLAST+** (`blastdbcmd`, `makeblastdb`) | all scripts | must be on `PATH`; tested with BLAST+ 2.14+ |
| **NCBI edirect** (`esearch`, `efetch`) | `make_subset_blastdb.*` | expected at `~/edirect`; install with `sh <(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)` |
| **NCBI taxonomy files** (`taxdb.btd`, `taxdb.bti`) | `subset_blastdb_by_accession.*` | must be present in the working directory or on `BLASTDB` path for taxid lookups to work |

**Note:** `BLAST Database error: Database memory map file error` (exit code 3) from
`blastdbcmd` is a misleading message that most commonly means the database path is wrong
or does not exist — check the path before suspecting a memory issue.

---

## Design notes

- **Standalone, not alias.** `blastdb_aliastool` alias databases require the full parent
  database to be present on every machine that uses them. These scripts produce
  self-contained databases that can be copied anywhere.
- **Accessions, not GI numbers.** NCBI no longer assigns GI numbers to new records.
  Accession.version strings are the stable, long-term identifier.
- **Nucleotide only.** Both scripts are hardcoded for `-dbtype nucl`.
