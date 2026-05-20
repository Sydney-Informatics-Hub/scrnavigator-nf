#!/usr/bin/env Rscript --vanilla
# Create 6 per-sample mouse Seurat test fixtures from E-MTAB-7324
# (Pijuan-Sala et al. 2019) via MouseGastrulationData::WTChimeraData().
#
# Sample mapping (WTChimeraData):
#   Samples 1, 3, 5  →  WT uninjected  (ctrl)
#   Samples 2, 4, 6  →  tdTomato       (stim)
#
# Outputs written to tests/data/rds/:
#   mouse_ctrl1.rds, mouse_ctrl2.rds, mouse_ctrl3.rds
#   mouse_stim1.rds, mouse_stim2.rds, mouse_stim3.rds
#
# REQUIRES: MouseGastrulationData, Seurat (>= 5), SingleCellExperiment, Matrix
# This script is NOT run inside the scrnavigator-nf-annotate container.
# Run once in a Bioconductor R session; commit the produced .rds files.
#
# Usage:
#   Rscript --vanilla tests/data/make_mouse_samples.R

library(Matrix)
library(SingleCellExperiment)
library(Seurat)

# TODO: Add BiocManager::install("MouseGastrulationData") if not installed
# Need to module load imagemagick on Gadi
library(MouseGastrulationData)

OUT_DIR          <- "tests/data/rds"
CELLS_PER_SAMPLE <- 500 # Gets you 1.9~2MB per sample
TOP_N_GENES      <- 2000
MIN_CELLS        <- 3

# WTChimeraData sample ID → (output file, condition)
sample_map <- list(
  "1" = list(file = "mouse_ctrl1.rds", condition = "ctrl"),
  "3" = list(file = "mouse_ctrl2.rds", condition = "ctrl"),
  "5" = list(file = "mouse_ctrl3.rds", condition = "ctrl"),
  "2" = list(file = "mouse_stim1.rds", condition = "stim"),
  "4" = list(file = "mouse_stim2.rds", condition = "stim"),
  "6" = list(file = "mouse_stim3.rds", condition = "stim")
)

expected <- file.path(OUT_DIR, vapply(sample_map, `[[`, character(1), "file"))
if (all(file.exists(expected))) {
  message("All mouse sample fixtures already exist, skipping.")
  quit(save = "no")
}

# ---------------------------------------------------------------------------
# 1. Download source data (ExperimentHub; cached after first run)
# ---------------------------------------------------------------------------
set.seed(42)
message("Downloading WTChimeraData samples 1-6 ...")
sce <- WTChimeraData(samples = 1:6)
message(sprintf("Loaded: %d genes x %d cells across %d samples",
                nrow(sce), ncol(sce), length(unique(sce$sample))))

# ---------------------------------------------------------------------------
# 2. Compute per-sample gene universe: top TOP_N_GENES by total UMI, >= MIN_CELLS
# ---------------------------------------------------------------------------
# Use the full dataset to rank genes so all samples share the same gene set
counts_full <- counts(sce)
gene_ncells  <- Matrix::rowSums(counts_full > 0)
gene_total   <- Matrix::rowSums(counts_full)
pass         <- which(gene_ncells >= MIN_CELLS)
top_genes    <- pass[order(gene_total[pass], decreasing = TRUE)][
  seq_len(min(TOP_N_GENES, length(pass)))
]
message(sprintf("Gene universe: %d genes (top %d by UMI, min %d cells)",
                length(top_genes), TOP_N_GENES, MIN_CELLS))

# ---------------------------------------------------------------------------
# 3. Per-sample: stratified subsample → filter genes → coerce → save
# ---------------------------------------------------------------------------
for (sid in names(sample_map)) {
  dest      <- file.path(OUT_DIR, sample_map[[sid]]$file)
  condition <- sample_map[[sid]]$condition
  sample_id <- tools::file_path_sans_ext(basename(dest))

  if (file.exists(dest)) {
    message("Already exists, skipping: ", dest)
    next
  }

  set.seed(42 + as.integer(sid))

  sce_s <- sce[, sce$sample == as.integer(sid)]
  message(sprintf("Sample %s (%s): %d cells available", sid, condition, ncol(sce_s)))

  # Stratified cell subsample by celltype.mapped
  meta      <- as.data.frame(colData(sce_s))
  strat_col <- if ("celltype.mapped" %in% colnames(meta)) "celltype.mapped" else NULL

  if (!is.null(strat_col) && ncol(sce_s) > CELLS_PER_SAMPLE) {
    ct_counts <- table(meta[[strat_col]])
    n_types   <- length(ct_counts)
    per_type  <- max(1L, floor(CELLS_PER_SAMPLE / n_types))
    keep_idx  <- unlist(lapply(names(ct_counts), function(ct) {
      idx <- which(meta[[strat_col]] == ct)
      sample(idx, min(per_type, length(idx)))
    }))
    keep_idx <- sample(keep_idx, min(CELLS_PER_SAMPLE, length(keep_idx)))
  } else {
    keep_idx <- sample(seq_len(ncol(sce_s)), min(CELLS_PER_SAMPLE, ncol(sce_s)))
  }

  sce_sub <- sce_s[top_genes, keep_idx]

  # Coerce to Seurat v5
  so <- as.Seurat(sce_sub, counts = "counts", data = NULL)
  so <- RenameAssays(so, originalexp = "RNA")
  DefaultAssay(so) <- "RNA"

  so$sample_id  <- sample_id
  so$condition  <- condition
  so$orig.ident <- sample_id

  saveRDS(so, dest)
  message(sprintf("Saved %s — %d cells x %d genes (%.1f MB)",
                  dest, ncol(so), nrow(so), file.info(dest)$size / 1e6))
}

# Sample 1 (ctrl): 2882 cells available
# Renaming default assay from originalexp to RNA
# Saved tests/data/rds/mouse_ctrl1.rds — 173 cells x 2000 genes (0.9 MB)
# Sample 3 (ctrl): 2468 cells available
# Renaming default assay from originalexp to RNA
# Saved tests/data/rds/mouse_ctrl2.rds — 161 cells x 2000 genes (0.8 MB)
# Sample 5 (ctrl): 2411 cells available
# Renaming default assay from originalexp to RNA
# Saved tests/data/rds/mouse_ctrl3.rds — 162 cells x 2000 genes (0.9 MB)
# Sample 2 (stim): 2597 cells available
# Renaming default assay from originalexp to RNA
# Saved tests/data/rds/mouse_stim1.rds — 180 cells x 2000 genes (1.0 MB)
# Sample 4 (stim): 1821 cells available
# Renaming default assay from originalexp to RNA
# Saved tests/data/rds/mouse_stim2.rds — 185 cells x 2000 genes (1.0 MB)
# Sample 6 (stim): 1047 cells available
# Renaming default assay from originalexp to RNA
# Saved tests/data/rds/mouse_stim3.rds — 170 cells x 2000 genes (0.9 MB)