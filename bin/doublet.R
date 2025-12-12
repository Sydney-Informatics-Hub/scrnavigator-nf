#!/usr/bin/env -S Rscript --vanilla
library(Seurat)
library(tidyverse)
library(DoubletFinder)

options(future.globals.maxSize = 1000*1024^2)

# Function definitions
find_min_pc <- function(seurat_obj) {
  stdvs <- seurat_obj@reductions$pca@stdev
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

  return(min(co1, co2))
}

run_doubletfinder_custom <- function(seurat_obj, cluster_res_name, multiplet_rate) {
  # From https://biostatsquid.com/doubletfinder-tutorial/
  min_pc <- find_min_pc(seurat_obj)

  # pK identification (no ground-truth)
  # introduces artificial doublets in varying props, merges with real data set and
  # preprocesses the data + calculates the prop of artficial neighrest neighbours,
  # provides a list of the proportion of artificial nearest neighbours for varying
  # combinations of the pN and pK
  sweep_list <- DoubletFinder::paramSweep(seurat_obj, PCs = 1:min_pc, sct = T)
  sweep_stats <- DoubletFinder::summarizeSweep(sweep_list)
  # computes a metric to find the optimal pK value (max mean variance normalised by modality coefficient)
  bcmvn <- DoubletFinder::find.pK(sweep_stats)

  # Optimal pK is the max of the bimodality coefficient (BCmvn) distribution
  optimal_pk <- bcmvn %>%
    dplyr::filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK)
  optimal_pk <- as.numeric(as.character(optimal_pk[[1]]))

  # Homotypic doublet % estimate
  annotations <- seurat_obj@meta.data[[cluster_res_name]]  # use the clusters as the user-defined cell types
  homotypic_prop <- DoubletFinder::modelHomotypic(annotations)  # get proportions of homotypic doublets
  nExp_poi <- round(multiplet_rate * nrow(seurat_obj@meta.data))  # multiply by number of cells to get the number of expected multiplets
  nExp_poi_adj <- round(nExp_poi * (1 - homotypic_prop))  # expected number of doublets

  # run DoubletFinder
  seurat_doublets <- DoubletFinder::doubletFinder(
    seu = seurat_obj,
    PCs = 1:min_pc,
    pN = 0.25,  # default
    pK = optimal_pk,  # the neighborhood size used to compute the number of artificial nearest neighbours
    nExp = nExp_poi_adj,  # number of expected real doublets
    reuse.pANN = NULL,
    sct = TRUE
  )

  # change name of metadata column with Singlet/Doublet information
  colnames(seurat_doublets@meta.data)[grepl('DF.classifications.*', colnames(seurat_doublets@meta.data))] <- "doublet_finder"

  return(seurat_doublets)
}

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

# Get number of cells and multiplet rate
# If no multiplet rate is defined, estimate it as 8e-6 per cell
n_cells <- nrow(so@meta.data)
multiplet_rate <- as.numeric(params$value[match("multiplet_rate", params$param)])
if (is.na(multiplet_rate)) {
  multiplet_rate <- 8e-6 * n_cells
}

# Detect doublets
cluster_res <- as.numeric(params$value[match("res", params$param)])
cluster_res_name <- paste0("SCT_snn_res.", cluster_res)
so_doublets_detected <- run_doubletfinder_custom(so, cluster_res_name, multiplet_rate)

# Remove doublets
so_doublets_removed <- subset(so_doublets_detected, subset = (doublet_finder == "Singlet"))

# Produce QC plots and summaries
p_doublet_umap <- DimPlot(so_doublets_detected, reduction = "umap", group.by = "doublet_finder")

doublet_summary_for_plotting <- so_doublets_detected@meta.data %>%
  dplyr::select(all_of(cluster_res_name), doublet_finder) %>%
  group_by(across(all_of(cluster_res_name)), doublet_finder) %>%
  summarise(n = n())

p_doublets_per_cluster <- doublet_summary_for_plotting %>%
  ggplot(aes(x = .data[[cluster_res_name]], y = n, fill = doublet_finder)) +
  geom_col() +
  theme_light()

p_doublet_props_per_cluster <- doublet_summary_for_plotting %>%
  pivot_wider(names_from = doublet_finder, values_from = n, values_fill = 0) %>%
  summarise(prop_doublets = Doublet / (Singlet + Doublet)) %>%
  ggplot(aes(x = .data[[cluster_res_name]], y = prop_doublets, fill = prop_doublets)) +
  geom_col() +
  theme_light() +
  ylim(0, 1) +
  scale_fill_continuous(type = "viridis", limits = c(0, 1)) +
  ylab("Proportion of doublets")

doublet_summary <- table(so_doublets_detected$doublet_finder)
doublet_summary_table <- tibble(
  sample = sample_id,
  n_singlets = doublet_summary[["Singlet"]],
  n_doublets = doublet_summary[["Doublet"]]
)

ggsave(paste0("qc_results/", sample_id, ".doublet_umap.", cluster_res, ".png"), p_doublet_umap)
ggsave(paste0("qc_results/", sample_id, ".doublets_per_cluster.", cluster_res, ".png"), p_doublets_per_cluster)
ggsave(paste0("qc_results/", sample_id, ".doublet_proportions_per_cluster.", cluster_res, ".png"), p_doublet_props_per_cluster)
doublet_summary_table %>% write_csv(paste0("qc_results/", sample_id, ".doublet_summary.csv"))

# Save Seurat object to file
SaveSeuratRds(so_doublets_detected, paste0(sample_id, ".doublets_detected.rds"))
SaveSeuratRds(so_doublets_removed, paste0(sample_id, ".doublets_removed.rds"))
