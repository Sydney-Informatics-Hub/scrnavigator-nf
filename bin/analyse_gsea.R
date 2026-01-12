#!/usr/bin/env -S Rscript --vanilla
library(Seurat)
library(tidyverse)
library(ggplot2)

# Get commandline arguments
args <- commandArgs(trailingOnly = TRUE)

cohort_id <- args[1]
all_rds_paths_file <- args[2:length(args)]

# Read in GSEA results from RDS files
full_gsea_results <- lapply(all_rds_paths_file, readRDS) %>% bind_rows()

# Recalculate FDR values across all tests
full_gsea_results <- full_gsea_results %>%
  mutate(FDR = p.adjust(pValue, method = "BH"))

# Reduce results down to one representative gene set per cluster
reduced_gsea_results <- full_gsea_results %>%
  filter(geneSet == top_gene_set)

# Plot bar graphs of results
reduced_gsea_results_split <- reduced_gsea_results %>%
  group_by(ref_group, test_group) %>%
  group_split()

gsea_plots <- lapply(reduced_gsea_results_split, function(x) {
  x %>%
    mutate(
      direction = case_when(normalizedEnrichmentScore > 0 ~ "Upregulated", .default = "Downregulated")
    ) %>%
    ggplot(aes(x = reorder(description, normalizedEnrichmentScore), y = normalizedEnrichmentScore, fill = direction)) +
    geom_bar(stat = "identity") +
    xlab(element_blank()) +
    ylab("log2 Enrichment Ratio") +
    theme(
      axis.title.x = element_text(size = 16),
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12)
    ) +
    geom_hline(yintercept = 0, linewidth = 1.5) +
    scale_fill_manual(name = "Direction", values = c(Upregulated = "orange", Downregulated = "royalblue2")) +
    coord_flip()
})
names(gsea_plots) <- sapply(reduced_gsea_results_split, function(x) {
  ref_group <- x$ref_group[1]
  test_group <- x$test_group[1]
  paste0(test_group, "_vs_", ref_group)
})

# Save plots to file
dir.create("plots")
for (n in names(gsea_plots)) {
  p <- gsea_plots[[n]]
  ggsave(paste0("plots/", cohort_id, ".gsea.", n, ".png"), p)
}

# Save results to file
write_csv(full_gsea_results, paste(cohort_id, "gsea.full.csv", sep = "."))
write_csv(reduced_gsea_results, paste(cohort_id, "gsea.full.reduced.csv", sep = "."))