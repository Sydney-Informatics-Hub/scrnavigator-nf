#!/usr/bin/env -S Rscript --vanilla
library(Seurat)
library(tidyverse)
library(ggplot2)

# Get commandline arguments
args <- commandArgs(trailingOnly = TRUE)

cohort_id <- args[1]
rds_path <- args[2]
grouping_fields <- args[3]

# Split grouping fields into vector
grouping_fields <- str_split_1(grouping_fields, ",")

# Read in Seurat object from RDS file
integrated <- readRDS(rds_path)

# Run pseudobulking
pseudo <- AggregateExpression(
  integrated,
  assays = "RNA",
  return.seurat = TRUE,
  group.by = c("orig.ident", grouping_fields)
)

# Set up comparison groupings
comparison_groups <- pseudo@meta.data %>%
  dplyr::select(all_of(grouping_fields)) %>%
  unite(comparison_group, sep = "_")
pseudo$comparison_group <- comparison_groups$comparison_group

Idents(pseudo) <- "comparison_group"

# Re-normalise the RNA count data
pseudo <- NormalizeData(pseudo, assay = "RNA", verbose = TRUE)

# Save pre-processed data to file
SaveSeuratRds(pseudo, paste0(cohort_id, ".pseudobulk.rds"))

# Save a table of comparison groups and the number of samples for each group
comparison_group_table <- table(pseudo$comparison_group) %>% as.data.frame()
colnames(comparison_group_table) <- c("comparison_group", "n_samples")
comparison_group_table %>%
  write_csv(paste0(cohort_id, ".comparison_groups.txt"))