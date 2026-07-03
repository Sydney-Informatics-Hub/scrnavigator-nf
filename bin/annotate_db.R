#!/usr/bin/env -S Rscript --vanilla
library(Seurat)
library(tidyverse)
library(SingleR)
library(ggplot2)
library(ggrepel)
library(scuttle)

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

# Perform annotation
min_cells_for_annotation <- as.integer(annotation_params$value[match("min_cells_for_annotation", annotation_params$param)])
annotation_db_file <- annotation_params$value[match("annotation_db", annotation_params$param)]
data_is_ensembl <- all(startsWith(rownames(integrated@assays$RNA), "ENS"))

use_hpca <- FALSE
use_mouse <- FALSE
if (!is.na(annotation_db_file)) {
  annotation_db <- readRDS(annotation_db_file)
} else if (species == 'human') {
  use_hpca <- TRUE
  if (data_is_ensembl) {
    annotation_db <- celldex::HumanPrimaryCellAtlasData(ensembl = TRUE)
  } else {
    annotation_db <- celldex::HumanPrimaryCellAtlasData(ensembl = FALSE)
  }
} else if (species == 'mouse') {
  use_mouse <- TRUE
  if (data_is_ensembl) {
    annotation_db <- celldex::MouseRNAseqData(ensembl = TRUE)
  } else {
    annotation_db <- celldex::MouseRNAseqData(ensembl = FALSE)
  }
} else {
  stop("Error: For species other than human and mouse, valid external annotation database RDS file must be provided.")
}

sce <- as.SingleCellExperiment(integrated, assay = "RNA")
sceM <- logNormCounts(sce)

if (use_hpca || use_mouse) {
  # Add the main-level annotations
  field_to_plot <- ifelse(use_hpca, "SingleR.hpca_main", "SingleR.mouse_main")
  predicted <- SingleR(test = sceM, ref = annotation_db, labels = annotation_db$label.main)
  keep <- table(predicted$labels) > min_cells_for_annotation
  integrated[[field_to_plot]] <- ifelse(keep[predicted$labels], predicted$labels, "Other")

  # Also add the fine-level annotations
  field_to_plot <- ifelse(use_hpca, "SingleR.hpca_fine", "SingleR.mouse_fine")
  predicted <- SingleR(test = sceM, ref = annotation_db, labels = annotation_db$label.fine)
  keep <- table(predicted$labels) > min_cells_for_annotation
  integrated[[field_to_plot]] <- ifelse(keep[predicted$labels], predicted$labels, "Other")
} else {
  predicted <- SingleR(test = sceM, ref= annotation_db, labels = annotation_db$label)
  keep <- table(predicted$labels) > min_cells_for_annotation
  integrated$SingleR.annotation <- ifelse(keep[predicted$labels], predicted$labels, "Other")
  field_to_plot <- "SingleR.annotation"
}

# Plot cell type assignment
default_res_name <- Misc(integrated, slot = "default_resolution")
p_singler_annotation_umap <- DimPlot(
  integrated,
  reduction = "umap",
  group.by = c(default_res_name, field_to_plot),
  label = TRUE,
  repel = TRUE,
  label.box = TRUE
)
ggsave(paste0("qc_results/", cohort_id, ".umap.singleR_annotation.png"), p_singler_annotation_umap)

# Plot cell type proportions per cluster
p_singler_annotation_prop <- integrated@meta.data %>%
  dplyr::select(any_of(c(default_res_name, field_to_plot))) %>%
  group_by(.data[[default_res_name]], .data[[field_to_plot]]) %>%
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
  ggtitle("Cluster cell type composition") +
  labs(x = "Cell Type", y = "Cluster", colour = "Prop. Cells", size = "Prop. Cells") +
  scale_colour_gradientn(breaks = c(0, 0.5, 1), colours = c("lightblue","beige","red"))
ggsave(paste0("qc_results/", cohort_id, ".cell_type_proportions.singleR_annotation.png"), p_singler_annotation_prop)

# Generate cel type proportion tables
if (use_hpca) {
  integrated@meta.data %>%
    dplyr::select(any_of(default_res_name), SingleR.hpca_main) %>%
    group_by(.data[[default_res_name]], SingleR.hpca_main) %>%
    summarise(n_cells = n()) %>%
    pivot_wider(names_from = default_res_name, values_from = n_cells) %>%
    write_csv(paste0("qc_results/", cohort_id, ".cell_type_proportions.singleR_annotation.hpca_main.csv"))

  integrated@meta.data %>%
    dplyr::select(any_of(default_res_name), SingleR.hpca_fine) %>%
    group_by(.data[[default_res_name]], SingleR.hpca_fine) %>%
    summarise(n_cells = n()) %>%
    pivot_wider(names_from = default_res_name, values_from = n_cells) %>%
    write_csv(paste0("qc_results/", cohort_id, ".cell_type_proportions.singleR_annotation.hpca_fine.csv"))

  available_annotations <- c("SingleR.hpca_main", "SingleR.hpca_fine")
} else if (use_mouse) {
  integrated@meta.data %>%
    dplyr::select(any_of(default_res_name), SingleR.mouse_main) %>%
    group_by(.data[[default_res_name]], SingleR.mouse_main) %>%
    summarise(n_cells = n()) %>%
    pivot_wider(names_from = default_res_name, values_from = n_cells) %>%
    write_csv(paste0("qc_results/", cohort_id, ".cell_type_proportions.singleR_annotation.mouse_main.csv"))

  integrated@meta.data %>%
    dplyr::select(any_of(default_res_name), SingleR.mouse_fine) %>%
    group_by(.data[[default_res_name]], SingleR.mouse_fine) %>%
    summarise(n_cells = n()) %>%
    pivot_wider(names_from = default_res_name, values_from = n_cells) %>%
    write_csv(paste0("qc_results/", cohort_id, ".cell_type_proportions.singleR_annotation.mouse_fine.csv"))

  available_annotations <- c("SingleR.mouse_main", "SingleR.mouse_fine")
} else {
  integrated@meta.data %>%
    dplyr::select(any_of(c(default_res_name, "SingleR.annotation"))) %>%
    group_by(.data[[default_res_name]], .data[["SingleR.annotation"]]) %>%
    summarise(n_cells = n()) %>%
    pivot_wider(names_from = default_res_name, values_from = n_cells) %>%
    write_csv(paste0("qc_results/", cohort_id, ".cell_type_proportions.singleR_annotation.csv"))

  available_annotations <- c("SingleR.annotation")
}

# Save pre-processed data to file
SaveSeuratRds(integrated, paste0(cohort_id, ".annotated.database.rds"))

# Save available annotations to file
sink("available_annotations.txt")
cat(paste(available_annotations, collapse = "\n"))
cat("\n")
sink()
