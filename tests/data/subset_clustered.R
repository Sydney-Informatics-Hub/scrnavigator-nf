#!/usr/bin/env Rscript --vanilla
# Apply clustering to the integrated test fixture.
# Run from scrnavigator-nf/:
#   singularity exec $SIF Rscript --vanilla tests/data/subset_clustered.R

here::i_am("tests/data/subset_clustered.R")

library(Seurat)

src  <- here::here("tests/data/rds/cohort.integrated.test.rds")
dest <- here::here("tests/data/rds/cohort.integrated.clustered.test.rds")

if (file.exists(dest)) {
  message("Fixture already exists: ", dest)
  quit(save = "no")
}

integrated <- readRDS(src)
message(sprintf("Loaded: %d cells, %d features", ncol(integrated), nrow(integrated)))

# Calculate min_pc (same logic as cluster_integrated.R)
vars <- integrated@reductions$pca@stdev
percent_var <- (vars / sum(vars)) * 100
cumulative <- cumsum(percent_var)
co1 <- which(cumulative > 90 & percent_var < 5)[1]
co2 <- sort(
  which(
    (percent_var[1:length(percent_var) - 1] - percent_var[2:length(percent_var)]) > 0.1
  ),
  decreasing = TRUE
)[1] + 1
min_pc <- min(co1, co2)

# Use pca (test fixture does not retain integrated.dr); set SCT as default assay
# so FindNeighbors names the graph SCT_snn, matching the column names annotate_custom.R expects
DefaultAssay(integrated) <- "SCT"
integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:min_pc, verbose = FALSE)
integrated <- FindClusters(integrated, resolution = c(0.5, 1), algorithm = 1, random.seed = 0, verbose = 0)

default_res_name <- "SCT_snn_res.1"
Idents(integrated) <- default_res_name
Misc(integrated, slot = "default_resolution") <- default_res_name

saveRDS(integrated, dest)
message(sprintf("Saved %s (%.1f MB)", dest, file.info(dest)$size / 1e6))
