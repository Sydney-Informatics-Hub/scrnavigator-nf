#!/usr/bin/env -S Rscript --vanilla
library(Seurat)
library(tidyverse)
library(ggplot2)

# Get commandline arguments
args <- commandArgs(trailingOnly = TRUE)

cohort_id <- args[1]
rds_path <- args[2]
ref_group <- args[3]
test_group <- args[4]

# Read in Seurat object from RDS file
pseudo <- readRDS(rds_path)

# Check that the ref and test groups are present
all_comparison_groups <- table(pseudo@meta.data$comparison_group)
stopifnot(ref_group %in% names(all_comparison_groups))
stopifnot(test_group %in% names(all_comparison_groups))

# Check that the ref and test groups have at least 2 samples each
stopifnot(all_comparison_groups[[ref_group]] >= 2)
stopifnot(all_comparison_groups[[test_group]] >= 2)
min_cells_per_group <- min(all_comparison_groups[[ref_group]], all_comparison_groups[[test_group]])

# Run differential expression
de <- FindMarkers(
  pseudo,
  ident.1 = test_group,
  ident.2 = ref_group,
  test.use = "DESeq2",
  min.cells.group = min_cells_per_group
)

# Add columns identifying the comparison
de$cohort <- cohort_id
de$test_group <- test_group
de$ref_group <- ref_group

# Save differential expression data to file
write_csv(de, paste(cohort_id, "de", test_group, ref_group, "csv", sep = "."))
saveRDS(de, paste(cohort_id, "de", test_group, ref_group, "Rds", sep = "."))
