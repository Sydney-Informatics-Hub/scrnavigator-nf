#!/usr/bin/env -S Rscript --vanilla
library(Seurat)
library(tidyverse)
library(ggplot2)

# Get commandline arguments
args <- commandArgs(trailingOnly = TRUE)

cohort_id <- args[1]
full_rds_paths_file <- args[2]
reduced_rds_paths_file <- args[3]

# Read in ORA results from RDS files
full_rds_paths <- scan(full_rds_paths_file, character())
reduced_rds_paths <- scan(reduced_rds_paths_file, character())
full_ora_results <- lapply(full_rds_paths, readRDS) %>% bind_rows()
reduced_ora_results <- lapply(reduced_rds_paths, readRDS) %>% bind_rows()

# TODO: Recalculate FDR values across all tests

# TODO: Plot bar graphs of results
# ora_plots <- lapply(names(ora_reduced_results), function(n) {
#   x <- ora_reduced_results[[n]]
#   x %>%
#     mutate(
#       log2EnrichmentRatio = log2(enrichmentRatio),
#       enrichment = case_when(log2EnrichmentRatio > 0 ~ "Enriched", .default = "Depleted")
#     ) %>%
#     ggplot(aes(x = reorder(description, log2EnrichmentRatio), y = log2EnrichmentRatio, fill = enrichment)) +
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
#     scale_fill_manual(name = "Enrichment", values = c(Enriched = "orange", Depleted = "royalblue2")) +
#     coord_flip()
# })
# names(ora_plots) <- names(ora_reduced_results)

# ora_plot_sizes <- lapply(names(ora_plots), function(n) {
#   p <- ora_plots[[n]]
#   df <- ora_reduced_results[[n]]
#   p_height <- 5 + dim(df)[1] * 0.2
#   p_width <- 5 + max(nchar(df$description)) * 0.1
#   return(list(height = p_height, width = p_width))
# })
# names(ora_plot_sizes) <- names(ora_plots)

# Save differential expression data to file
write_csv(full_ora_results, paste(cohort_id, "ora.full.csv", sep = "."))
write_csv(reduced_ora_results, paste(cohort_id, "ora.full.reduced.csv", sep = "."))
