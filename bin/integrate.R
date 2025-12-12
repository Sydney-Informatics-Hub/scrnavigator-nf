#!/usr/bin/env -S Rscript --vanilla
library(Seurat)
library(tidyverse)
library(ggplot2)
library(clustree)

options(future.globals.maxSize = 1000*1024^2)

# Get commandline arguments
args <- commandArgs(trailingOnly = TRUE)

cohort_id <- args[1]
rds_paths <- args[2:length(args)]

# Read in Seurat objects from RDS files
all_so <- lapply(rds_paths, readRDS)

# Create QC output directory
dir.create("qc_results")

# TODO: Merge, integrate, cluster
integrated <- all_so[[1]]  # Placeholder

# Save Seurat object to file
SaveSeuratRds(integrated, paste0(cohort_id, ".integrated.rds"))
