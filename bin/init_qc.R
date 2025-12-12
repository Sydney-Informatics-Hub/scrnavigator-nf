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

# Plot nCount_RNA, nFeature_RNA, and percent.mt distributions
sample_metadata <- so@meta.data %>%
  rownames_to_column(var = "barcode")

dist_cols <- c("nCount_RNA", "nFeature_RNA")
if ("percent.mt" %in% colnames(sample_metadata)) {
  dist_cols <- c(dist_cols, "percent.mt")
}

p_dist <- sample_metadata %>%
  pivot_longer(cols = dist_cols) %>%
  ggplot(aes(x = orig.ident, y = value, fill = orig.ident, colour = orig.ident)) +
  facet_wrap(~ name, scales = "free_y") +
  theme_light() +
  theme(legend.position = "none")

p_dist_log <- p_dist + scale_y_log10()

p_dist <- p_dist + geom_violin(alpha = 0.1)
p_dist_log <- p_dist_log + geom_violin(alpha = 0.1)

ggsave(paste0("qc_results/", sample_id, ".count_distributions.png"), p_dist)
ggsave(paste0("qc_results/", sample_id, ".count_distributions.log.png"), p_dist_log)

# Generate pre-filtering nFeature_RNA vs nCount_RNA plot
p_fvc_pre_filtered <- sample_metadata %>%
  ggplot(aes(x = nCount_RNA, y = nFeature_RNA, colour = percent.mt)) +
  geom_point(size = 0.1) +
  scale_color_viridis_c() +
  scale_x_log10() +
  scale_y_log10() +
  theme_light() +
  annotation_logticks(colour = "lightgrey")

# Apply hard filtering thresholds
min_ncount <- as.numeric(params$value[match("min_ncount", params$param)])
max_ncount <- as.numeric(params$value[match("max_ncount", params$param)])
min_nfeature <- as.numeric(params$value[match("min_nfeature", params$param)])
max_nfeature <- as.numeric(params$value[match("max_nfeature", params$param)])
min_mt_pct <- as.numeric(params$value[match("min_mt_pct", params$param)])
max_mt_pct <- as.numeric(params$value[match("max_mt_pct", params$param)])

# Handle NAs and validate
if (is.na(min_ncount)) {
  min_ncount <- 0
}
stopifnot(min_ncount >= 0)

if (is.na(max_ncount)) {
  max_ncount <- Inf
}
stopifnot(max_ncount >= 0)

if (is.na(min_nfeature)) {
  min_nfeature <- 0
}
stopifnot(min_nfeature >= 0)

if (is.na(max_nfeature)) {
  max_nfeature <- Inf
}
stopifnot(max_nfeature >= 0)

if (is.na(min_mt_pct)) {
  min_mt_pct <- 0
}
stopifnot(min_mt_pct >= 0 && min_mt_pct <= 100)

if (is.na(max_mt_pct)) {
  max_mt_pct <- 100
}
stopifnot(max_mt_pct >= 0 && max_mt_pct <= 100)

so <- so
so@meta.data <- so@meta.data %>%
  mutate(
    qc_threshold = case_when(
      nCount_RNA >= min_ncount &
        nCount_RNA <= max_ncount &
        nFeature_RNA >= min_nfeature &
        nFeature_RNA <= max_nfeature &
        percent.mt >= min_mt_pct &
        percent.mt <= max_mt_pct ~ "keep",
      .default = "remove"
    )
  )
so_filtered <- subset(so, subset = (qc_threshold == "keep"))

# Remove any additional cells that user specifies within file provided to the cells_to_remove parameter
cells_to_remove_file <- params$value[match("cells_to_remove", params$param)]
if (!is.na(cells_to_remove_file) && cells_to_remove_file != "") {
  cells_to_remove <- scan(cells_to_remove_file, character())
  rm_cells <- rownames(so_filtered@meta.data) %in% cells_to_remove
  keep_cells <- !rm_cells
  cells_to_keep <- rownames(so_filtered@meta.data[keep_cells,])
  so_filtered <- subset(so_filtered, cells = cells_to_keep)
}

# Generate post-filtering nFeature_RNA vs nCount_RNA plot
sample_metadata_filtered <- so_filtered@meta.data %>%
  rownames_to_column(var = "barcode")
p_fvc_post_filtered <- sample_metadata_filtered %>%
  ggplot(aes(x = nCount_RNA, y = nFeature_RNA, colour = percent.mt)) +
  geom_point(size = 0.1) +
  scale_color_viridis_c() +
  scale_x_log10() +
  scale_y_log10() +
  theme_light() +
  annotation_logticks(colour = "lightgrey")

ggsave(paste0("qc_results/", sample_id, ".feature_count_plot.pre_filtered.png"), p_fvc_pre_filtered)
ggsave(paste0("qc_results/", sample_id, ".feature_count_plot.post_filtered.png"), p_fvc_post_filtered)

# Generate a filter summary table
n_cells_pre_filtered <- nrow(so@meta.data)
n_cells_post_filtered <- nrow(so_filtered@meta.data)

filter_summary <- tibble(
  sample = sample_id,
  n_cells_pre_filtered = n_cells_pre_filtered,
  n_cells_post_filtered = n_cells_post_filtered
) %>%
  mutate(
    prop_filtered = 1 - n_cells_post_filtered / n_cells_pre_filtered
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

filter_summary %>% write_csv(paste0("qc_results/", sample_id, ".hard_filter_summary.csv"))

# Perform SCTransform and clustering
cluster_method <- params$value[match("cluster_method", params$param)]
stopifnot(cluster_method %in% c("louvain", "leiden"))
cluster_algorithm <- ifelse(cluster_method == "louvain", 1, 4)

cluster_resolutions <- params$value[match("resolutions", params$param)]
cluster_resolutions <- str_split_1(cluster_resolutions, ",") %>%
  sapply(as.numeric, USE.NAMES = FALSE)

random_seed <- ifelse(cluster_algorithm == 1, 0, 1)

# Run initial SCTransform and PCA
sct_filtered <-
  Seurat::SCTransform(
    so_filtered,
    vars.to.regress = c("percent.mt"),
    verbose = FALSE
  ) |>
  Seurat::RunPCA()

# Calculate the minimum number of PCs that explain the majority of the variation
stdvs <- sct_filtered@reductions$pca@stdev
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

sct_filtered <-
  Seurat::RunUMAP(sct_filtered, dims = 1:min_pc, verbose = FALSE) |>
  Seurat::FindNeighbors(dims = 1:min_pc, verbose = FALSE) |>
  FindClusters(
    resolution = cluster_resolutions,
    algorithm = cluster_algorithm,
    random.seed = random_seed,
    verbose = 0
  )

clustree_plot <- clustree::clustree(sct_filtered, prefix = "SCT_snn_res.") +
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
  p_fvc_clustered <- sct_filtered@meta.data %>%
    ggplot(aes(x = nCount_RNA, y = nFeature_RNA, col = percent.mt)) +
      geom_point(size = 0.3) +
      facet_wrap(res_named) +
      scale_x_log10() +
      scale_y_log10() +
      theme_light() +
      viridis::scale_color_viridis() +
      annotation_logticks(side = "lb", colour = "lightgrey") +
      ggtitle(paste0(sample_id, ": ", res_named))
  ggsave(paste0("qc_results/", sample_id, ".feature_count_plot.cluster.", res,".png"))

  dir.create(paste0("qc_results/cluster_cells/", res), recursive = TRUE)
  all_clusters <- unique(sct_filtered@meta.data[[res_named]])
  for (cls in all_clusters) {
    cluster_file <- paste0("qc_results/cluster_cells/", res, "/", cls, ".txt")
    cluster_cells <- rownames(sct_filtered@meta.data[sct_filtered@meta.data[[res_named]] == cls,])
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
Idents(sct_filtered) <- default_res_name

# Save Seurat object to file
SaveSeuratRds(sct_filtered, paste0(sample_id, ".qc.rds"))
