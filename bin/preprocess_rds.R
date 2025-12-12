#!/usr/bin/env -S Rscript --vanilla
library(Seurat)
library(tidyverse)

# Get commandline arguments
args <- commandArgs(trailingOnly = TRUE)

sample_id <- args[1]
rds_path <- args[2]
metadata_file <- args[3]

# Read in Seurat object from RDS file
so <- readRDS(rds_path)

# Read in sample metadata
metadata <- read_csv(metadata_file)

# Update metadata
so@project.name <- sample_id
so$orig.ident <- as.factor(sample_id)

# Add in additional sample metadata
for (field in metadata$field) {
  so[[field]] <- metadata$value[match(field, metadata$field)]
}

# Ensure both gene symbols and Ensembl IDs are present in the RNA assay metadata
# TODO: Support other species (at least mouse)
ens_rownames <- all(startsWith(rownames(so@assays$RNA), "ENS"))
if (ens_rownames) {
  ens_ids <- rownames(so@assays$RNA)
  gene_symbols <- AnnotationDbi::mapIds(
    EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
    keys = ens_ids,
    column = "SYMBOL",
    keytype = "GENEID"
  )
  stopifnot(all(ens_ids == names(gene_symbols)))
} else {
  gene_symbols <- rownames(so@assays$RNA)
  ens_ids <- AnnotationDbi::mapIds(
    EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
    keys = gene_symbols,
    column = "GENEID",
    keytype = "SYMBOL"
  )
  stopifnot(all(gene_symbols == names(ens_ids)))
}
have_ens_ids <- "gene_versions" %in% colnames(so@assays$RNA@meta.data) &&
  all(startsWith(as.character(so@assays$RNA@meta.data$gene_versions), "ENS"))
have_gene_symbols <- "gene_symbols" %in% colnames(so@assays$RNA@meta.data) &&
  !all(startsWith(as.character(so@assays$RNA@meta.data$gene_symbols), "ENS"))
if (!have_ens_ids) {
  so@assays$RNA@meta.data$gene_versions <- ens_ids
}
if (!have_gene_symbols) {
  so@assays$RNA@meta.data$gene_symbols <- gene_symbols
}

# Add mitochondrial gene percentage per cell
# TODO: Support other species (at least mouse)
mt_pattern <- "^MT-"

rn <- rownames(so@assays$RNA)
if (all(startsWith(rn, "ENS"))) {
  rownames(so@assays$RNA) <- so@assays$RNA@meta.data$gene_symbols
  so$percent.mt <- Seurat::PercentageFeatureSet(so, pattern = mt_pattern)
  rownames(so@assays$RNA) <- rn
} else {
  so$percent.mt <- Seurat::PercentageFeatureSet(so, pattern = mt_pattern)
}

# Save pre-processed data to file
SaveSeuratRds(so, paste0(sample_id, ".preprocessed.rds"))
