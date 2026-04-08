#!/usr/bin/env -S Rscript --vanilla
# https://jorainer.github.io/ensembldb/articles/ensembldb.html#getting-or-building-ensdb-databasespackages
library(ensembldb) # Create sqlitedb of EnsDb class and query it for annotation
library(AnnotationHub) # pre-built ensdb objects
library(AnnotationFilter)

# Get commandline arguments
args <- commandArgs(trailingOnly = TRUE)

species <- args[1] # e.g. "Homo sapiens", "Mus musculus"
ref_version <- args[2] # ensembl version e.g. latest v113
db_cache <- args[3] # Where temporary and final db files are saved

# For ref_version, either
# Select the annotation version that matches the version of transcriptome
# See: https://www.10xgenomics.com/support/software/cell-ranger/latest/release-notes/cr-reference-release-notes#2024-a 
# or 
# the latest EnsDb version

# Downloads few ~MB cache on first run. This is the list of entries across species
ah <- AnnotationHub(cache = db_cache)

# Query for human EnsDb objects
ens_query <- query(ah, c("Ensembl", ref_version, species))
ens_query

# AnnotationHub with 1 record
# snapshotDate(): 2025-10-29
# names(): AH113665
# $dataprovider: Ensembl
# $species: Homo sapiens
# $rdataclass: EnsDb
# $rdatadateadded: 2023-04-25
# $title: Ensembl 110 EnsDb for Homo sapiens
# $description: Gene and protein annotations for Homo sapiens based on Ensembl version 110.
# $taxonomyid: 9606
# $genome: GRCh38
# $sourcetype: ensembl
# $sourceurl: http://www.ensembl.org
# $sourcesize: NA
# $tags: c("110", "Annotation", "AnnotationHubSoftware", "Coverage", "DataImport", "EnsDb", "Ensembl",
#   "Gene", "Protein", "Sequencing", "Transcript") 
# retrieve record with 'object[["AH113665"]]' 

stopifnot("Multiple db entries found. Cannot automatically select." = length(ens_query) == 1)
db_uid <- names(ens_query@.db_uid) # target db release that matches ref transcriptome v.

# Downloads db to cache. Takes 20-30 mins on first run.
ensdb <- ah[[db_uid]]
ensdb

# EnsDb for Ensembl:
# |Backend: SQLite
# |Db type: EnsDb
# |Type of Gene ID: Ensembl Gene ID
# |Supporting package: ensembldb
# |Db created by: ensembldb package from Bioconductor
# |script_version: 0.3.10
# |Creation time: Mon Aug  7 09:02:07 2023
# |ensembl_version: 110
# |ensembl_host: 127.0.0.1
# |Organism: Homo sapiens
# |taxonomy_id: 9606
# |genome_build: GRCh38
# |DBSCHEMAVERSION: 2.2
# |common_name: human
# |species: homo_sapiens
# | No. of genes: 71440.
# | No. of transcripts: 278545.
# |Protein data available.
# Export to standalone .sqlite

cached_sql_path <- dbconn(ensdb)@dbname
# Unreadable file name when downloaded:
#".../.cache/AnnotationHub/f43f73d14706_120411"
# Move and rename to "./EnsDb_homo_sapiens_v110.sqlite"
species_idx <- which(metadata(ensdb)[["name"]] == "species")
species_name <- metadata(ensdb)[species_idx, 2]

clean_db_path <- paste0(
  "./EnsDb_",
  species_name,
  "_v",
  ensemblVersion(ensdb),
  ".sqlite"
) 

# Move/rename cached sqlite file
file.rename(from = cached_sql_path, to = clean_db_path)
file.info(clean_db_path)$size / 1e6  # Verify, size in MB

# Load saved ensdb
#test_load <- EnsDb(clean_db_path)
