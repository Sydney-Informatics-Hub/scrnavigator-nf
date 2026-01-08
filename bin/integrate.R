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

# Get sample names
sample_names <- sapply(
  all_so,
  function(so) {
    so$orig.ident[1]
  },
  USE.NAMES = FALSE
)
names(all_so) <- sample_names

# Ensure all gene IDs are the same (symbols or Ensembl IDs)
using_ens_ids <- sapply(all_so, function(so) {
  all(startsWith(rownames(so@assays$RNA), "ENS"))
})

if (any(using_ens_ids) && !all(using_ens_ids)) {
  all_so <- lapply(all_so, function(so) {
    if (all(startsWith(rownames(so@assays$RNA), "ENS"))) {
      return(so)
    }
    if (any(duplicated(so@assays$RNA@meta.data$gene_versions))) {
      warning("Duplicated Ensembl IDs have been detected when converting from gene symbols. For each duplicated ID, the first will be kept and the rest will be dropped.")
      keep <- !duplicated(so@assays$RNA@meta.data$gene_versions)
      s <- s[keep, ]
    }
    if (any(is.na(so@assays$RNA@meta.data$gene_versions))) {
      warning("Missing Ensembl IDs have been detected when converting from gene symbols. These will be dropped.")
      keep <- !is.na(so@assays$RNA@meta.data$gene_versions)
      s <- s[keep, ]
    }
    rownames(so@assays$RNA) <- so@assays$RNA@meta.data$gene_versions
    return(so)
  })
}

# Merge data
first_so <- all_so[[1]]
remaining_so <- all_so[2:length(all_so)]

merged <- merge(
  first_so,
  remaining_so,
  add.cell.ids = sample_names
)

# Strip back to RNA assay
DefaultAssay(merged) <- "RNA"
merged <- DietSeurat(merged, assays = c("RNA"))

# Remove old cluster annotations
meta_cols <- colnames(merged@meta.data)
sct_cluster_cols <- startsWith(meta_cols, "SCT_snn_res.")
pann_cols <- startsWith(meta_cols, "pANN_")  # NOTE: Do we want to remove these columns?
seurat_cluster_col <- meta_cols == "seurat_clusters"
meta_cols_to_remove <- sct_cluster_cols | pann_cols | seurat_cluster_col
merged@meta.data <- merged@meta.data[, !meta_cols_to_remove]

# Perform SCTransform on merged data
# Run initial SCTransform and PCA
merged <-
  Seurat::SCTransform(
    merged,
    vars.to.regress = c("percent.mt"),
    verbose = FALSE
  ) |>
  Seurat::RunPCA()

# Calculate the minimum number of PCs that explain the majority of the variation
stdvs <- merged@reductions$pca@stdev
percent_stdv <- (stdvs / sum(stdvs)) * 100
cumulative <- cumsum(percent_stdv)
co1 <- which(cumulative > 90 & percent_stdv < 5)[1]
co2 <- sort(
  which(
    (
      percent_stdv[1:length(percent_stdv) - 1] -
        percent_stdv[2:length(percent_stdv)]
    ) > 0.1
  ),
  decreasing = TRUE
)[1] + 1
min_pc <- min(co1, co2)

merged <- Seurat::RunUMAP(merged, dims = 1:min_pc, verbose = FALSE)

# Generate merged but pre-integrated UMAP
p_pre_integrated_umap <- DimPlot(merged, reduction = "umap", group.by = c("orig.ident"))
ggsave(paste0("qc_results/", cohort_id, ".umap.merged.png"), p_pre_integrated_umap)

integrated <- IntegrateLayers(merged, method = CCAIntegration, normalization.method = "SCT")
integrated <- Seurat::RunUMAP(integrated, dims = 1:min_pc, reduction = "integrated.dr", verbose = FALSE)

# Generate integrated UMAP
p_post_integrated_umap <- DimPlot(integrated, reduction = "umap", group.by = c("orig.ident"))
ggsave(paste0("qc_results/", cohort_id, ".umap.integrated.png"), p_post_integrated_umap)

# Save Seurat object to file
SaveSeuratRds(integrated, paste0(cohort_id, ".integrated.rds"))

# Save gene symbols and Ensembl IDs to an RDS file
gene_symbols <- integrated@assays$RNA@meta.data %>%
  dplyr::select(gene_symbols, gene_versions)
saveRDS(gene_symbols, "gene_symbols.Rds")
