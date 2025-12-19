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
sig_de <- de %>%
  dplyr::filter(
    test_group == test_group_val,
    ref_group == ref_group_val,
    sig == "sig"
  )

# Prepare data for ORA
sig_de <- sig_de %>%
  mutate(abs_avg_log2FC = abs(avg_log2FC)) %>%
  arrange(desc(avg_log2FC))

background_genes <- sig_de %>%
  pull(Gene)

using_ens_ids <- all(startsWith(sig_de$Gene, "ENS"))
if (using_ens_ids) {
  gene_type <- "ensembl_gene_id"
} else {
  gene_type <- "genesymbol"
}

comparison <- paste0(test_group_val, "_vs_", ref_group_val)
dir.create(comparison, recursive = TRUE)

databases <- listGeneSet()

# Run ORA
WebGestaltR(
  enrichMethod = "ORA",
  interestGene = sig_de$Gene,
  interestGeneType = gene_type,
  referenceGene = background_genes,
  referenceGeneType = gene_type,
  enrichDatabase = databases$name[startsWith(databases$name, "pathway")],
  isOutput = TRUE,
  nThreads = 1,
  outputDirectory = comparison,
  projectName = comparison
)

# Get results files (if they exist)
results_file <- paste0(comparison, "/Project_", comparison, "/enrichment_results_", comparison, ".txt")
ap_file <- paste0(comparison, "/Project_", comparison, "/enriched_geneset_ap_clusters_", comparison, ".txt")
if (file.exists(results_file) && file.exists(ap_file)) {
  # Read results file
  ora_results <- read_tsv(results_file, show_col_types = FALSE)

  # Read clusters file
  ap_file <- paste0(comparison, "/Project_", comparison, "/enriched_geneset_ap_clusters_", comparison, ".txt")
  ap_raw <- readLines(ap_file) %>% str_split("\t")
  ap_clusters <- lapply(1:length(ap_raw), function(x) {
    n <- paste0("cluster_", x)
    data.frame(cluster = n, gene_sets = ap_raw[[x]])
  }) %>%
    do.call(rbind, .)

  # Add cluster information to results file
  top_ap_clusters <- ap_clusters %>%
    group_by(cluster) %>%
    dplyr::filter(row_number() == 1) %>%
    rename(top_gene_set = gene_sets)
  all_ap_clusters <- ap_clusters %>%
    group_by(cluster) %>%
    summarise(gene_sets_combined = paste(gene_sets, collapse = ";"))
  ap_clusters <- ap_clusters %>%
    left_join(all_ap_clusters, by = "cluster") %>%
    left_join(top_ap_clusters, by = "cluster")

  ora_results <- ora_results %>%
    left_join(ap_clusters, by = join_by(geneSet == gene_sets))

  # Add comparison details to dataframes
  ora_results$cohort <- cohort_id
  ora_results$test_group <- test_group
  ora_results$ref_group <- ref_group

  # Save to file
  write_csv(ora_results, paste(cohort_id, "ora", test_group, ref_group, "csv", sep = "."))
  saveRDS(ora_results, paste(cohort_id, "ora", test_group, ref_group, "Rds", sep = "."))
}
