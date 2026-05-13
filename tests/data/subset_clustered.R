#!/usr/bin/env Rscript --vanilla
# Apply clustering to the integrated test fixture.
# Run from scrnavigator-nf/:
#   singularity exec $SIF Rscript --vanilla tests/data/subset_clustered.R

library(Seurat)

src  <- "tests/data/rds/cohort.integrated.test.rds"
dest <- "tests/data/rds/cohort.integrated.clustered.test.rds"

if (file.exists(dest)) {
  message("Fixture already exists: ", dest)
  quit(save = "no")
}

integrated <- readRDS(src)
message(sprintf("Loaded: %d cells, %d features", ncol(integrated), nrow(integrated)))

# Calculate min_pc (same logic as cluster_integrated.R)
stdvs <- integrated@reductions$pca@stdev
percent_stdv <- (stdvs / sum(stdvs)) * 100
cumulative <- cumsum(percent_stdv)
co1 <- which(cumulative > 90 & percent_stdv < 5)[1]
co2 <- sort(
  which(
    (percent_stdv[1:length(percent_stdv) - 1] - percent_stdv[2:length(percent_stdv)]) > 0.1
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
