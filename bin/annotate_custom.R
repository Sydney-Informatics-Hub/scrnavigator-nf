#!/usr/bin/env -S Rscript --vanilla
library(Seurat)
library(tidyverse)
library(SingleR)
library(ggplot2)
library(ggrepel)

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

# Get species information
species <- annotation_params$value[match("species", annotation_params$param)]

# Read in custom marker gene sets
custom_marker_genes_file <- annotation_params$value[match("custom_marker_genes", annotation_params$param)]
if (!is.na(custom_marker_genes_file)) {
  custom_marker_genes <- read_csv(custom_marker_genes_file)
  custom_programs <- custom_marker_genes$cell_type
  custom_marker_genes <- strsplit(custom_marker_genes$gene_ids, ";")
  names(custom_marker_genes) <- custom_programs

  if(!(length(custom_marker_genes) > 1)) {
    stop("Error: You have only defined one custom gene program. Custom cell type annotation requires at least two gene programs. Please define additional cell types.")
  }
} else {
  stop("Error: No custom marker gene file provided.")
}

# Convert gene IDs to/from Ensembl if necessary
marker_genes_are_ensembl <- sapply(custom_marker_genes, function(s) {
  all(startsWith(s, "ENS"))
})
if (any(marker_genes_are_ensembl)) {
  stopifnot(all(marker_genes_are_ensembl))
  all_marker_genes_are_ensembl <- TRUE
} else {
  all_marker_genes_are_ensembl <- FALSE
}

ensdb_file <- annotation_params$value[match("ens_db", annotation_params$param)]

if (!is.na(ensdb_file)) {
  endsb <- readRDS(ensdb_file)
} else if (species == "human") {
  ensdb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
} else if (species == "mouse") {
  ensdb <- EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79
} else {
  stop("Error: Species is neither human nor mouse and no EnsDb database was provided.")
}

data_is_ensembl <- all(startsWith(rownames(integrated@assays$RNA), "ENS"))

if (data_is_ensembl && !all_marker_genes_are_ensembl) {
  custom_marker_genes <- lapply(custom_marker_genes, function(g) {
    gene2ens <- AnnotationDbi::mapIds(
      ensdb,
      keys = g,
      column = "GENEID",
      keytype = "SYMBOL"
    )
    keep <- !is.na(gene2ens)
    return(gene2ens[keep])
  })
} else if (!data_is_ensembl && all_marker_genes_are_ensembl) {
  custom_marker_genes <- lapply(custom_marker_genes, function(g) {
    ens2gene <- AnnotationDbi::mapIds(
      ensdb,
      keys = g,
      column = "SYMBOL",
      keytype = "GENEID"
    )
    keep <- !is.na(ens2gene)
    return(ens2gene[keep])
  })
}

# Add gene set module scores
integrated <- AddModuleScore(
  integrated,
  features = custom_marker_genes,
  name = names(custom_marker_genes)
)

# Remove numeric suffix from cluster names
numeric_suffix_clusters <- paste0(names(custom_marker_genes), 1:length(custom_marker_genes))
stopifnot(all(numeric_suffix_clusters %in% colnames(integrated@meta.data)))  # Sanity check
colnames(integrated@meta.data)[colnames(integrated@meta.data) %in% numeric_suffix_clusters] <- names(custom_marker_genes)
stopifnot(all(names(custom_marker_genes) %in% colnames(integrated@meta.data)))  # Sanity check

# Summarise scores for each gene program
default_res_name <- Misc(integrated, slot = "default_resolution")
custom_scores <- integrated@meta.data %>%
  dplyr::select(Cluster = all_of(default_res_name), all_of(custom_programs)) %>%
  pivot_longer(
    cols = all_of(custom_programs),
    names_to = "cell_type",
    values_to = "score"
  ) %>%
  group_by(Cluster, cell_type) %>%
  summarise(
    avg_score    = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    sd_score     = sd(score, na.rm = TRUE),
    n_cells      = n(),
    .groups      = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))
write_csv(custom_scores, paste0("qc_results/", cohort_id, ".custom_cell_type_scores.csv"))

# Plot scores for each gene program
p_custom_annotations_cluster_heatmap <- ggplot(custom_scores, aes(x = Cluster, y = cell_type, fill = median_score)) +
  geom_tile() +
  theme_light() +
  ggtitle("Cluster Scores") +
  scale_fill_gradientn(colours = c("lightblue","beige","red"))
ggsave(paste0("qc_results/", cohort_id, ".heatmap.custom_annotations_per_cluster.png"), p_custom_annotations_cluster_heatmap)

# Predict the cell type using the maximum score
# For 2 or 3 cell types, the maximum score is used
# For 4+ cell types, the maximum score is used IF it is significantly greater than the next highest score
# Significance is determined using the median absolute difference (MAD) of all scores for a cell and
# a MAD threshold supplied by the user (default is 1.0)
custom_annotation_mad_threshold <- as.numeric(annotation_params$value[match("custom_annotation_mad_threshold", annotation_params$param)])

argmax_scores <- apply(integrated@meta.data[custom_programs], 1, which.max) %>% unlist
integrated$custom_cell_type.max_score <- custom_programs[argmax_scores]

# Also calculate how much the top two scores differ
# Can only calculate when more than 3 cell programs are present
# Determines how many MADs the top 2 scores differ by
# If they differ by less than 1 MAD, mark as ambiguous
if (length(custom_programs) > 3) {
  s <- integrated@meta.data[custom_programs]
  integrated$custom_cell_type.top_2_score_mad_diff <- apply(s, 1, function(x) {
    xv <- unlist(x)
    mad_x <- stats::mad(xv, constant = 1)
    top_2_x <- xv[order(xv, decreasing = TRUE)][1:2]
    return((top_2_x[1] - top_2_x[2]) / mad_x)
  }) %>% unlist
  
  # If the top two scores are too close, mark the cell as ambiguous
  integrated@meta.data <- integrated@meta.data %>% mutate(
    custom_cell_type = case_when(
      custom_cell_type.top_2_score_mad_diff >= custom_annotation_mad_threshold ~ custom_cell_type.max_score,
      .default = "Ambiguous"
    )
  )

  available_annotations <- c("custom_cell_type.max_score", "custom_cell_type")
} else {
  integrated@meta.data <- integrated@meta.data %>% mutate(
    custom_cell_type = custom_cell_type.max_score
  )

  available_annotations <- c("custom_cell_type")
}

# Plot the cell type assignments on the UMAP
p_custom_annotation_umap <- DimPlot(
  integrated,
  reduction = "umap",
  group.by = "custom_cell_type",
  label = TRUE,
  repel = TRUE,
  label.box = TRUE
)
ggsave(paste0("qc_results/", cohort_id, ".umap.custom_annotation.png"), p_custom_annotation_umap)

# Plot assignments per cluster
p_custom_annotations_per_cluster <- integrated@meta.data %>%
  dplyr::select(any_of(default_res_name), custom_cell_type) %>%
  group_by(.data[[default_res_name]], custom_cell_type) %>%
  summarise(n_cells = n()) %>%
  group_by(.data[[default_res_name]]) %>%
  mutate(total_cells = sum(n_cells)) %>%
  purrr::set_names(c("cluster", "cell_type", "n_cells", "total_cells")) %>%
  mutate(prop_cells = n_cells / total_cells) %>%
  ggplot(aes(x = cell_type, y = cluster, colour = prop_cells, size = prop_cells)) +
    geom_point() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(
        fill = NA,
        colour = "grey70",
        linewidth = rel(1)
    ),
    panel.grid = element_line(colour = "grey87"),
    panel.grid.major = element_line(linewidth = rel(0.5)),
    panel.grid.minor = element_line(linewidth = rel(0.25)),
    axis.ticks = element_line(colour = "grey70", linewidth = rel(0.5)),
    strip.background = element_rect(
      fill = "grey70",
      colour = NA
    ),
    complete = TRUE
  ) +
  ggtitle("Cluster cell type composition (including ambiguous cell types)") +
  labs(x = "Cell Type", y = "Cluster", colour = "Prop. Cells", size = "Prop. Cells") +
  scale_colour_gradientn(breaks = c(0, 0.5, 1), colours = c("lightblue","beige","red"))
ggsave(paste0("qc_results/", cohort_id, ".cell_type_proportions.custom_annotation.png"), p_custom_annotations_per_cluster)

integrated@meta.data %>%
  dplyr::select(any_of(default_res_name), custom_cell_type) %>%
  group_by(.data[[default_res_name]], custom_cell_type) %>%
  summarise(n_cells = n()) %>%
  pivot_wider(names_from = default_res_name, values_from = n_cells) %>%
  write_csv(paste0("qc_results/", cohort_id, ".cell_type_proportions.custom_annotation.csv"))

# Plot the module scores per program
for (ct in custom_programs) {
  p_custom_annotation_scores <- FeaturePlot(
    integrated,
    reduction = "umap",
    features = c(ct)
  ) +
    scale_colour_gradient2(
      low = 'royalblue',
      mid = 'lightgrey',
      high = 'indianred2',
      midpoint = 0
    )
  ggsave(paste0("qc_results/", cohort_id, ".custom_annotation_scores.", ct,".png"), p_custom_annotation_scores)
}

# Save pre-processed data to file
SaveSeuratRds(integrated, paste0(cohort_id, ".annotated.custom.rds"))

# Save available annotations to file
sink("available_annotations.txt")
cat(paste(available_annotations, collapse = "\n"))
cat("\n")
sink()
