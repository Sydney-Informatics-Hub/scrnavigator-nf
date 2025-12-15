#!/usr/bin/env -S Rscript --vanilla
library(Seurat)
library(tidyverse)

# Get commandline arguments
args <- commandArgs(trailingOnly = TRUE)

sample_id <- args[1]
rds_path <- args[2]
species_file <- args[3]

# Read in Seurat object from RDS file
integrated <- readRDS(rds_path)

# Read in sample metadata
metadata <- read_csv(metadata_file)

# Read in species parameters
species_params <- read_csv(species_file)

# Create QC output directory
dir.create("qc_results")

# Get species information
species <- species_params$value[match("species", species_params$param)]

# TODO: Perform annotation

# Save pre-processed data to file
SaveSeuratRds(integrated, paste0(sample_id, ".annotated.database.rds"))
