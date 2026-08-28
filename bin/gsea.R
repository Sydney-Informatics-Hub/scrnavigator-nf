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
species <- args[5]
opt_db_file <- args[6]

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

# Get organism ID for ORA
organism_lookup <- c(
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
if (species %in% names(organism_lookup)) {
  organism <- organism_lookup[species]
} else if (organism %in% listOrganism()) {
  organism <- species
} else {
  organism <- "others"
}

# Get database for ORA
if (!is.na(opt_db_file) && opt_db_file != "") {
  enrichment_database <- NULL
  enrichment_database_file <- opt_db_file
} else {
  databases <- listGeneSet()
  enrichment_database <- databases$name[grepl("^geneontology_.*_noRedundant$", databases$name, perl = TRUE)]
  enrichment_database_file <- NULL
}

# Run GSEA
WebGestaltR(
  enrichMethod = "GSEA",
  organism = organism,
  interestGene = comparison_de,
  interestGeneType = gene_type,
  referenceGene = background_genes,
  referenceGeneType = gene_type,
  enrichDatabase = enrichment_database,
  enrichDatabaseFile = enrichment_database_file,
  isOutput = TRUE,
  nThreads = 1,
  outputDirectory = comparison,
  projectName = "gsea"
)

# Get results files (if they exist)
results_file <- paste0(comparison, "/Project_gsea/enrichment_results_gsea.txt")
ap_file <- paste0(comparison, "/Project_gsea/enriched_geneset_ap_clusters_gsea.txt")
if (file.exists(results_file) && file.exists(ap_file)) {
  # Read results file
  gsea_results <- read_tsv(results_file, show_col_types = FALSE)

  # Read clustesr file
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

  gsea_results <- gsea_results %>%
    left_join(ap_clusters, by = join_by(geneSet == gene_sets))

  # Add comparison details to dataframes
  gsea_results$cohort <- cohort_id
  gsea_results$test_group <- test_group_val
  gsea_results$ref_group <- ref_group_val

  # Save to file
  write_csv(gsea_results, paste(cohort_id, "gsea", test_group_val, ref_group_val, "csv", sep = "."))
  saveRDS(gsea_results, paste(cohort_id, "gsea", test_group_val, ref_group_val, "Rds", sep = "."))
}
