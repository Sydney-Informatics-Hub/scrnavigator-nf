#!/usr/bin/env -S Rscript --vanilla
library(Seurat)
library(tidyverse)
library(ggplot2)
library(clustree)

options(future.globals.maxSize = 1000*1024^2)

# Get commandline arguments
args <- commandArgs(trailingOnly = TRUE)

cohort_id <- args[1]
rds_path <- args[2]
params_file <- args[3]

# Read in Seurat object from RDS file
integrated <- readRDS(rds_path)

# Read in parameters
params <- read_csv(params_file)

# Create QC output directory
dir.create("qc_results")

# Perform clustering
cluster_method <- params$value[match("cluster_method", params$param)]
stopifnot(cluster_method %in% c("louvain", "leiden"))
cluster_algorithm <- ifelse(cluster_method == "louvain", 1, 4)

cluster_resolutions <- params$value[match("resolutions", params$param)]
cluster_resolutions <- str_split_1(cluster_resolutions, ",") %>%
  sapply(as.numeric, USE.NAMES = FALSE)

random_seed <- ifelse(cluster_algorithm == 1, 0, 1)

# Calculate the minimum number of PCs that explain the majority of the variation
vars <- integrated@reductions$pca@stdev ^ 2
percent_var <- (vars / sum(vars)) * 100
cumulative <- cumsum(percent_var)
co1 <- which(cumulative > 90 & percent_var < 5)[1]
co2 <- sort(
  which(
    (
      percent_var[1:length(percent_var) - 1] -
        percent_var[2:length(percent_var)]
    ) > 0.1
  ),
  decreasing = TRUE
)[1] + 1
min_pc <- min(co1, co2)

integrated <- FindNeighbors(integrated, reduction = "integrated.dr", dims = 1:min_pc, verbose = FALSE)
integrated <- FindClusters(
  integrated,
  resolution = cluster_resolutions,
  algorithm = cluster_algorithm,
  random.seed = random_seed,
  verbose = 0
)

# Additionally correct the SCT counts after integration
integrated <- PrepSCTFindMarkers(integrated)

# Run clustree
clustree_plot <- clustree::clustree(integrated, prefix = "SCT_snn_res.") +
  ggtitle(cohort_id) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.justification.right = "top"
  )
ggsave(paste0("qc_results/", cohort_id, ".clustree.png"), clustree_plot)

# Set a default cluster resolution
default_res <- 1
default_res_param <- as.numeric(params$value[match("integrated_resolution", params$param)])
if (!is.na(default_res_param)) {
  default_res <- default_res_param
}
default_res_name <- paste0("SCT_snn_res.", default_res)
Idents(integrated) <- default_res_name
Misc(integrated, slot = "default_resolution") <- default_res_name

# Plot UMAP with clusters at integrated resolution
p_integrated_clustered_umap <- DimPlot(integrated, reduction = "umap", group.by = c("orig.ident", default_res_name), label = TRUE)
ggsave(paste0("qc_results/", cohort_id, ".umap.integrated.clusters.png"), p_integrated_clustered_umap)

# Plot per-cluster nFeature_RNA vs nCount_RNA plot at integrated resolution
p_fvc_integrated_clustered <- integrated@meta.data %>%
  ggplot(aes(x = nCount_RNA, y = nFeature_RNA, col = percent.mt)) +
    geom_point(size = 0.3) +
    facet_wrap(default_res_name) +
    scale_x_log10() +
    scale_y_log10() +
    theme_light() +
    viridis::scale_color_viridis() +
    annotation_logticks(side = "lb", colour = "lightgrey") +
    ggtitle(paste0("Integrated Dataset: ", default_res_name))
ggsave(paste0("qc_results/", cohort_id, ".feature_count_plot.integrated.cluster.png"), p_fvc_integrated_clustered)

# Get number of clusters
n_clusters <- table(integrated@meta.data[default_res_name])
n_clusters <- tibble(data.frame(n_clusters))
colnames(n_clusters) <- c("cluster", "n_cells")
n_clusters %>% write_csv(paste0("qc_results/", cohort_id, ".integrated.clusters.csv"))

# Save Seurat object to file
SaveSeuratRds(integrated, paste0(cohort_id, ".integrated.clustered.rds"))
