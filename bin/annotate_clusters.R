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
manual_annotations_file <- args[4]  # May be NA

# Read in Seurat object from RDS file
integrated <- readRDS(rds_path)

# Read in sample metadata
metadata <- read_csv(metadata_file)

# Read in annotation parameters
annotation_params <- read_csv(annotation_file)

# Create QC output directory
dir.create("qc_results")

# Get species information
species <- annotation_params$value[match("species", annotation_params$param)]

if (!is.na(manual_annotations_file)) {
  cluster_cell_type_assignments <- read_csv(manual_annotations_file)
} else {
  # Get user-supplied annotation to use or set default
  cluster_annotation <- annotation_params$value[match("cluster_annotation", annotation_params$param)]
  if (is.na(cluster_annotation)) {
    annotation_priority <- c(
      "custom_cell_type",
      "custom_cell_type.max_score",
      "SingleR.annotation",
      "SingleR.hpca_main",
      "SingleR.hpca_fine",
      "Phase"
    )
    cluster_annotation <- which(annotation_priority %in% colnames(integrated@meta.data))[1]
  } else {
    stopifnot(cluster_annotation %in% colnames(integrated@meta.data))
  }

  # Get cell type proportion threshold
  cell_type_proportion_threshold <- as.numeric(annotation_params$value[match("cell_type_proportion_threshold", annotation_params$param)])

  # Call the cell type of the cluster
  default_res_name <- Idents(integrated)
  cluster_annotations <- integrated@meta.data[c(default_res_name, cluster_annotation)]
  colnames(cluster_annotations) <- c("cluster", "cell_type")
  cluster_annotations$cell_id <- rownames(cluster_annotations)
  cluster_sizes <- cluster_annotations %>%
    as_tibble %>%
    group_by(cluster) %>%
    summarise(cluster_size = n())
  cluster_annotations_summary <- cluster_annotations %>%
    as_tibble %>%
    group_by(cluster, cell_type) %>%
    summarise(n_cells = n()) %>%
    left_join(cluster_sizes, by = "cluster") %>%
    mutate(prop_cells = n_cells / cluster_size)

  # Create a dot plot of the cell type annotations per cluster
  p_cluster_cell_type_proportions <- cluster_annotations_summary %>%
    arrange(prop_cells) %>%
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
    ggtitle("Cluster cell type composition") +
    labs(x = "Cell Type", y = "Cluster", colour = "Prop. Cells", size = "Prop. Cells") +
    scale_colour_gradientn(breaks = c(0, 0.5, 1), colours = c("lightblue","beige","red"))
  ggsave(paste0("qc_results/", cohort_id, ".cluster_cell_type_proportions.png"), p_cluster_cell_type_proportions)

  # Write the cell type annotations per cluster to a CSV file
  cluster_annotations_summary %>%
    mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
    write_csv(paste0("qc_results/", cohort_id, ".cluster_cell_type_proportions.csv"))

  # Determine the consensus annotations
  cluster_cell_type_assignments <- cluster_annotations_summary %>%
    group_by(cluster) %>%
    dplyr::filter(prop_cells == max(prop_cells)) %>%
    mutate(
      cell_type_consensus = case_when(
        prop_cells >= cell_type_proportion_threshold ~ cell_type,
        .default = "Ambiguous"
      )
    ) %>%
    dplyr::select(cluster, cell_type_consensus) %>%
    ungroup() %>%
    dplyr::rename(cell_type = cell_type_consensus) %>%
    write_csv(paste0(cohort_id, ".cluster_cell_type_assignments.csv"))
}

# Assign the cell type annotations to each cluster
integrated@meta.data$cluster_annotation <- cluster_cell_type_assignments$cell_type[match(integrated@meta.data[[default_res_name]], cluster_cell_type_assignments$cluster)]

# Plot cluster assignments
p_cluster_cell_type_assignments_umap <- DimPlot(
  integrated,
  reduction = "umap",
  group.by = "cluster_annotation",
  label = TRUE,
  repel = TRUE,
  label.box = TRUE
)
ggsave(paste0("qc_results/", cohort_id, ".umap.cluster_cell_type_assignments.png"), p_cluster_cell_type_assignments_umap)

# Save pre-processed data to file
SaveSeuratRds(integrated, paste0(cohort_id, ".annotated.clusters.rds"))
