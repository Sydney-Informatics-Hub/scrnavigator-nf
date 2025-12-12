#!/usr/bin/env -S Rscript --vanilla
library(Seurat)
library(tidyverse)
library(ggplot2)
library(clustree)

options(future.globals.maxSize = 1000*1024^2)

# Get commandline arguments
args <- commandArgs(trailingOnly = TRUE)

sample_id <- args[1]
rds_path <- args[2]
params_file <- args[3]

# Read in Seurat object from RDS file
so <- readRDS(rds_path)

# Read in parameters
params <- read_csv(params_file)

# Create QC output directory
dir.create("qc_results")

# Perform SCTransform and clustering
cluster_method <- params$value[match("cluster_method", params$param)]
stopifnot(cluster_method %in% c("louvain", "leiden"))
cluster_algorithm <- ifelse(cluster_method == "louvain", 1, 4)

cluster_resolutions <- params$value[match("resolutions", params$param)]
cluster_resolutions <- str_split_1(cluster_resolutions, ",") %>%
  sapply(as.numeric, USE.NAMES = FALSE)

random_seed <- ifelse(cluster_algorithm == 1, 0, 1)

# Run initial SCTransform and PCA
sct <-
  Seurat::SCTransform(
    so,
    vars.to.regress = c("percent.mt"),
    verbose = FALSE
  ) |>
  Seurat::RunPCA()

# Calculate the minimum number of PCs that explain the majority of the variation
stdvs <- sct@reductions$pca@stdev
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

sct <-
  Seurat::RunUMAP(sct, dims = 1:min_pc, verbose = FALSE) |>
  Seurat::FindNeighbors(dims = 1:min_pc, verbose = FALSE) |>
  FindClusters(
    resolution = cluster_resolutions,
    algorithm = cluster_algorithm,
    random.seed = random_seed,
    verbose = 0
  )

clustree_plot <- clustree::clustree(sct, prefix = "SCT_snn_res.") +
  ggtitle(sample_id) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.justification.right = "top"
  )

ggsave(paste0("qc_results/", sample_id, ".clustree.png"), clustree_plot)

# Plot per-cluster nFeature_RNA vs nCount_RNA plots for each cluster resolution
# Also print out lists of cell/barcode IDs for each cluster for each resolution
for (res in cluster_resolutions) {
  res_named <- paste0("SCT_snn_res.", res)
  p_fvc_clustered <- sct@meta.data %>%
    ggplot(aes(x = nCount_RNA, y = nFeature_RNA, col = percent.mt)) +
      geom_point(size = 0.3) +
      facet_wrap(res_named) +
      scale_x_log10() +
      scale_y_log10() +
      theme_light() +
      viridis::scale_color_viridis() +
      annotation_logticks(side = "lb", colour = "lightgrey") +
      ggtitle(paste0(sample_id, ": ", res_named))
  ggsave(paste0("qc_results/", sample_id, ".feature_count_plot.cluster.", res,".png"), p_fvc_clustered)

  dir.create(paste0("qc_results/cluster_cells/", res), recursive = TRUE)
  all_clusters <- unique(sct@meta.data[[res_named]])
  for (cls in all_clusters) {
    cluster_file <- paste0("qc_results/cluster_cells/", res, "/", cls, ".txt")
    cluster_cells <- rownames(sct@meta.data[sct@meta.data[[res_named]] == cls,])
    sink(cluster_file)
    cat(paste(cluster_cells, collapse = "\n"))
    sink()
  }
}

# Set a default cluster resolution
default_res <- 1
default_res_param <- as.numeric(params$value[match("res", params$param)])
if (!is.na(default_res_param)) {
  default_res <- default_res_param
}
default_res_name <- paste0("SCT_snn_res.", default_res)
Idents(sct) <- default_res_name

# Save Seurat object to file
SaveSeuratRds(sct, paste0(sample_id, ".sct_clustered.rds"))
