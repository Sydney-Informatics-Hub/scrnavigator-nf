#!/usr/bin/env -S Rscript --vanilla
library(Seurat)
library(tidyverse)
library(ggplot2)
library(WebGestaltR)

# Get commandline arguments
args <- commandArgs(trailingOnly = TRUE)

cohort_id <- args[1]
rds_path <- args[2]
ref_group_val <- args[3]
test_group_val <- args[4]

# Read in DE results
de <- readRDS(rds_path)

# Filter the current comparison and for significatn DEGs
comparison_de <- de %>%
  dplyr::filter(
    test_group == test_group_val,
    ref_group == ref_group_val
  )

# Prepare data for GSEA
comparison_de <- comparison_de %>%
  arrange(desc(avg_log2FC)) %>%
  dplyr::select(Gene, avg_log2FC)

background_genes <- comparison_de %>%
  pull(Gene)

using_ens_ids <- all(startsWith(comparison_de$Gene, "ENS"))
if (using_ens_ids) {
  gene_type <- "ensembl_gene_id"
} else {
  gene_type <- "genesymbol"
}

comparison <- paste0(test_group_val, "_vs_", ref_group_val)
dir.create(comparison, recursive = TRUE)

databases <- listGeneSet()

# Run GSEA
WebGestaltR(
  enrichMethod = "GSEA",
  interestGene = comparison_de,
  interestGeneType = gene_type,
  referenceGene = background_genes,
  referenceGeneType = gene_type,
  enrichDatabase = databases$name[grepl("^geneontology_.*_noRedundant$", databases$name, perl = TRUE)],
  isOutput = TRUE,
  nThreads = 1,
  outputDirectory = comparison,
  projectName = comparison
)

# Get results files (if they exist)
results_file <- paste0(comparison, "/Project_", comparison, "/enrichment_results_", comparison, ".txt")
if (file.exists(results_file)) {
  gsea_results <- read_tsv(results_file, show_col_types = FALSE)

  # Reduce to top gene set per cluster
  ap_file <- paste0(comparison, "/Project_", comparison, "/enriched_geneset_ap_clusters_", comparison, ".txt")
  ap_raw <- readLines(ap_file) %>% str_split("\t")
  ap_clusters <- lapply(1:length(ap_raw), function(x) {
    n <- paste0("cluster_", x)
    data.frame(cluster = n, gene_sets = ap_raw[[x]])
  }) %>%
    do.call(rbind, .) %>%
    group_by(cluster) %>%
    dplyr::filter(row_number() == 1)

  gsea_reduced_results <- gsea_results %>%
    left_join(ap_clusters, by = join_by(geneSet == gene_sets)) %>%
    drop_na(cluster)

  # Add comparison details to dataframes
  gsea_results$cohort <- cohort_id
  gsea_results$test_group <- test_group
  gsea_results$ref_group <- ref_group
  gsea_reduced_results$cohort <- cohort_id
  gsea_reduced_results$test_group <- test_group
  gsea_reduced_results$ref_group <- ref_group

  # Save to file
  write_csv(gsea_results, paste(cohort_id, "gsea", test_group, ref_group, "csv", sep = "."))
  saveRDS(gsea_results, paste(cohort_id, "gsea", test_group, ref_group, "Rds", sep = "."))
  write_csv(gsea_reduced_results, paste(cohort_id, "gsea", test_group, ref_group, "reduced.csv", sep = "."))
  saveRDS(gsea_reduced_results, paste(cohort_id, "gsea", test_group, ref_group, "reduced.Rds", sep = "."))
}
