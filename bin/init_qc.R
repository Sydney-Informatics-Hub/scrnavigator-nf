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
min_ncount <- params$value[match("min_ncount", params$param)]
max_ncount <- params$value[match("max_ncount", params$param)]
min_nfeature <- params$value[match("min_nfeature", params$param)]
max_nfeature <- params$value[match("max_nfeature", params$param)]
min_mt_pct <- params$value[match("min_mt_pct", params$param)]
max_mt_pct <- params$value[match("max_mt_pct", params$param)]

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

# TODO: Remove cells that user supplies with --remove_cells parameter

# Perform SCTransform and clustering
cluster_method <- params$value[match("cluster_method", params$param)]
stopifnot(cluster_method %in% c("louvain", "leiden"))
cluster_algorithm <- ifelse(cluster_method == "louvain", 1, 4)

cluster_resolutions <- params$value[match("resolutions", params$param)]
cluster_resolutions <- str_split_1(cluster_resolutions) %>%
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
}

# TODO:
#   - Print out lists of all cell/barcode IDs for each cluster into separate files
#   - Allow the user to supply a parameter --remove_cells file1[,file2[...]]
#   - The user can look at the QC data generated by this script and provide the files
#     for the clusters they want to remove, as well as any other cells they might want to remove

# TODO: Set a default cluster resolution