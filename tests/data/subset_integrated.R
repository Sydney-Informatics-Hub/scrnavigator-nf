#!/usr/bin/env Rscript --vanilla
# Subset the integrated RDS to a small test fixture (~10 MB).
# Run from scrnavigator-nf/:
#   singularity exec exec $SINGULARITY_CACHEDIR/sydneyinformaticshub-scrnavigator-nf-annotate.img Rscript --vanilla tests/data/subset_integrated.R

library(Seurat)

# cohort.integrated.rds is produced by running the pipeline through integration
# on the test RDS files. See docs/how-to-run-tests.md — "Cell cycle RDS" section.
src  <- "tests/data/rds/cohort.integrated.rds"
dest <- "tests/data/rds/cohort.integrated.test.rds"

if (file.exists(dest)) {
  message("Fixture already exists: ", dest)
  quit(save = "no")
}

integrated <- readRDS(src)
message(sprintf("Loaded: %d cells, %d features", ncol(integrated), nrow(integrated)))

set.seed(42)
cells_keep <- sample(colnames(integrated), min(500, ncol(integrated)))
subset_obj <- integrated[, cells_keep]

subset_obj <- DietSeurat(
  subset_obj,
  assays    = c("RNA", "SCT"),
  dimreducs = c("pca", "umap"),
  layers    = c("counts", "data")
)

saveRDS(subset_obj, dest)
message(sprintf("Saved %s (%.1f MB)", dest, file.info(dest)$size / 1e6))
