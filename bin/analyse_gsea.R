#!/usr/bin/env -S Rscript --vanilla
library(Seurat)
library(tidyverse)
library(ggplot2)

# Get commandline arguments
args <- commandArgs(trailingOnly = TRUE)

cohort_id <- args[1]
full_rds_paths_file <- args[2]
reduced_rds_paths_file <- args[3]

# Read in GSEA results from RDS files
full_rds_paths <- scan(full_rds_paths_file, character())
reduced_rds_paths <- scan(reduced_rds_paths_file, character())
full_gsea_results <- lapply(full_rds_paths, readRDS) %>% bind_rows()
reduced_gsea_results <- lapply(reduced_rds_paths, readRDS) %>% bind_rows()

# TODO: Recalculate FDR values across all tests

# TODO: Plot bar graphs of results
# gsea_plots <- lapply(names(gsea_reduced_results), function(n) {
#   x <- gsea_reduced_results[[n]]
#   x %>%
#     mutate(
#       direction = case_when(normalizedEnrichmentScore > 0 ~ "Upregulated", .default = "Downregulated")
#     ) %>%
#     ggplot(aes(x = reorder(description, normalizedEnrichmentScore), y = normalizedEnrichmentScore, fill = direction)) +
#     geom_bar(stat = "identity") +
#     xlab(element_blank()) +
#     ylab("log2 Enrichment Ratio") +
#     theme(
#       axis.title.x = element_text(size = 16),
#       axis.text.y = element_text(size = 12),
#       axis.text.x = element_text(size = 12),
#       legend.title = element_text(size = 12),
#       legend.text = element_text(size = 12)
#     ) +
#     geom_hline(yintercept = 0, linewidth = 1.5) +
#     scale_fill_manual(name = "Direction", values = c(Upregulated = "orange", Downregulated = "royalblue2")) +
#     coord_flip()
# })
# names(gsea_plots) <- names(gsea_reduced_results)

# gsea_plot_sizes <- lapply(names(gsea_plots), function(n) {
#   p <- gsea_plots[[n]]
#   df <- gsea_reduced_results[[n]]
#   p_height <- 5 + dim(df)[1] * 0.2
#   p_width <- 5 + max(nchar(df$description)) * 0.1
#   return(list(height = p_height, width = p_width))
# })
# names(gsea_plot_sizes) <- names(gsea_plots)

# Save differential expression data to file
write_csv(full_gsea_results, paste(cohort_id, "gsea.full.csv", sep = "."))
write_csv(reduced_gsea_results, paste(cohort_id, "gsea.full.reduced.csv", sep = "."))