#!/usr/bin/env -S Rscript --vanilla
# https://jorainer.github.io/ensembldb/articles/ensembldb.html#getting-or-building-ensdb-databasespackages
library(ensembldb, quietly = TRUE) # Create sqlitedb of EnsDb class and query it for annotation
library(AnnotationHub, quietly = TRUE) # pre-built ensdb objects

# Get commandline arguments
args <- commandArgs(trailingOnly = TRUE)

common_name <- args[1] # e.g. "human", "mouse"
ref_version <- args[2] # e.g. latest ensembl version "v113"
db_cache <- args[3] # Where temporary and final db files are saved ".cache"

# For ref_version, either
# Select the annotation version that matches the version of transcriptome
# See: https://www.10xgenomics.com/support/software/cell-ranger/latest/release-notes/cr-reference-release-notes#2024-a 
# or 
# the latest EnsDb version

# From bin/gsea.R. Get corresponding abbreviated name
abbrev_lookup <- c(
  human = "hsapiens",
  mouse = "mmusculus",
  rat = "rnorvegicus",
  chicken = "ggallus",
  dog = "cfamiliaris",
  cattle = "btaurus",
  cow = "btaurus",
  zebrafish = "drerio",
  pig = "sscrofa",
  fruit_fly = "dmelanogaster"
)

# Full species name spaced (e.g. for ensdb lookup)
species_lookup <- c(
  human = "homo sapiens",
  mouse = "mus musculus",
  rat = "rattus norvegicus",
  chicken = "gallus gallus",
  dog = "canis lupus_familiaris",
  cattle = "bos taurus",
  cow = "bos taurus",
  zebrafish = "danio rerio",
  pig = "sus scrofa",
  fruit_fly = "drosophila_melanogaster"
)

stopifnot("Species not found in lookup table." = common_name %in% names(abbrev_lookup))
abbreviated_name <- abbrev_lookup[common_name]
species_name <- species_lookup[common_name]

# Downloads few ~MB cache on first run. This is the list of entries across species
ah <- AnnotationHub(cache = db_cache)

# Query for human EnsDb objects
ens_query <- query(ah, c("Ensembl", ref_version, species_name))
ens_query

if (length(ens_query) != 1) {
  stopifnot("No db entries found. Check your common name and version inputs." = length(ens_query) < 1)
  stopifnot("Multiple db entries found. Cannot automatically select." = length(ens_query) > 1)
}
db_uid <- names(ens_query@.db_uid) # target db release that matches ref transcriptome v.

# Downloads db to cache. Takes 1-10 mins on first run.
ensdb <- ah[[db_uid]]
ensdb

# Unreadable file name when downloaded:
# .cache/AnnotationHub/f43f73d14706_120411"
# Move and rename to "./EnsDb_hsapiens_v110.sqlite"
cached_sql_path <- dbconn(ensdb)@dbname

clean_db_path <- paste0(
  "./EnsDb_",
  abbreviated_name,
  "_v",
  ensemblVersion(ensdb),
  ".sqlite"
) 

# Move/rename cached sqlite file
file.rename(from = cached_sql_path, to = clean_db_path)
file.info(clean_db_path)$size / 1e6  # Verify, size in MB

# Load saved ensdb
#test_load <- EnsDb(clean_db_path)
