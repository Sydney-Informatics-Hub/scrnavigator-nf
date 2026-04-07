#!/usr/bin/env -S Rscript --vanilla
library(Seurat)
library(tidyverse)
library(SingleR)
#BiocManager::install(c("ensembldb", "AnnotationHub", "AnnotationFilter"))
library(ensembldb) # Create sqlitedb of EnsDb class and query it for annotation
library(AnnotationHub) # pre-built ensdb objects
library(AnnotationFilter)

# Inputs
species <- "Homo sapiens"

# Select the annotation version that matches the version of transcriptome
# See: https://www.10xgenomics.com/support/software/cell-ranger/latest/release-notes/cr-reference-release-notes#2024-a 
ref_version <- c("Ensembl", "110")

# Downloads few ~MB cache on first run
ah <- AnnotationHub(cache = ".cache/AnnotationHub")

# Query for human EnsDb objects
ens_query <- query(ah, c(ref_version, species))
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

# Downloads db release into memory. Takes 20-30 mins
ensdb <- ah[[db_uid]]
ensdb


# Export to standalone .sqlite
db_path <- paste0("EnsDb_SPECIES_vENSVERSION.sqlite") #"EnsDb_Hsapiens_v113.sqlite"
ensdb_local <- copyDb(ensdb, destination = "./db_path")

# Verify saved on disk
file.info(db_path)$size / 1e6  # Size in MB

# Load saved ensdb
ensdb2 <- EnsDb(db_path)
ensdb2