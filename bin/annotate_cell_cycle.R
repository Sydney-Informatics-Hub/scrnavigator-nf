#!/usr/bin/env -S Rscript --vanilla
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ensembldb)

# Get commandline arguments
args <- commandArgs(trailingOnly = TRUE)

cohort_id <- args[1]
rds_path <- args[2]
annotation_file <- args[3]

# Read in Seurat object from RDS file
integrated <- readRDS(rds_path)

# Read in annotation parameters
annotation_params <- read_csv(annotation_file)

# Create QC output directory
dir.create("qc_results")

# Get cell cycle genes for species
species <- annotation_params$value[match("species", annotation_params$param)]
s_genes_file <- annotation_params$value[match("s2_genes", annotation_params$param)]
g2m_genes_file <- annotation_params$value[match("g2m_genes", annotation_params$param)]

if (file.exists(s_genes_file) && file.exists(g2m_genes_file)) {
  s_genes <- scan(s_genes_file, character())
  g2m_genes <- scan(g2m_genes_file, character())
} else if (species == "human") {
  s_genes <- cc.genes$s.genes
  g2m_genes <- cc.genes$g2m.genes
} else {
  stop("Error: For non-human species, valid S and G2M gene files must be provided.")
}

# Convert cell cycle genes between symbols and Ensembl IDs where necessary
ensdb_file <- annotation_params$value[match("ens_db", annotation_params$param)]

if (file.exists(ensdb_file)) {
  ensdb <- EnsDb(ensdb_file)
} else if (species == "human") {
  ensdb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
} else if (species == "mouse") {
  ensdb <- EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79
} else {
  stop("Error: Species is neither human nor mouse and no EnsDb database was provided.")
}

data_is_ensembl <- all(startsWith(rownames(integrated@assays$RNA), "ENS"))
s_genes_are_ensembl <- all(startsWith(s_genes, "ENS"))
g2m_genes_are_ensembl <- all(startsWith(g2m_genes, "ENS"))

if (data_is_ensembl) {
  if (!s_genes_are_ensembl) {
    s_genes <- AnnotationDbi::mapIds(
      ensdb,
      keys = s_genes,
      column = "GENEID",
      keytype = "SYMBOL"
    )
  }
  if (!g2m_genes_are_ensembl) {
    g2m_genes <- AnnotationDbi::mapIds(
      ensdb,
      keys = g2m_genes,
      column = "GENEID",
      keytype = "SYMBOL"
    )
  }
} else {
  if (s_genes_are_ensembl) {
    s_genes <- AnnotationDbi::mapIds(
      ensdb,
      keys = s_genes,
      column = "SYMBOL",
      keytype = "GENEID"
    )
  }
  if (g2m_genes_are_ensembl) {
    g2m_genes <- AnnotationDbi::mapIds(
      ensdb,
      keys = g2m_genes,
      column = "SYMBOL",
      keytype = "GENEID"
    )
  }
}

# Annotate cell cycle
integrated <- CellCycleScoring(
  integrated,
  s.features = s_genes,
  g2m.features = g2m_genes
)

# Plot cell cycle annotations
p_cell_cycle_umap <- DimPlot(integrated, reduction = "umap", group.by = "Phase")
ggsave(paste0("qc_results/", cohort_id, ".umap.cell_cycle.png"), p_cell_cycle_umap)

# Save pre-processed data to file
SaveSeuratRds(integrated, paste0(cohort_id, ".annotated.cell_cycle.rds"))

# Save available annotations to file
sink("available_annotations.txt")
cat("Phase\n")
sink()
