#!/usr/bin/env -S Rscript --vanilla
library(Seurat)
library(tidyverse)
library(ggplot2)

# Get commandline arguments
args <- commandArgs(trailingOnly = TRUE)

cohort_id <- args[1]
all_rds_paths <- args[2:length(args)]

# Read in ORA results from RDS files
full_ora_results <- lapply(all_rds_paths, readRDS) %>% bind_rows()

# Recalculate FDR values across all tests
full_ora_results <- full_ora_results %>%
  mutate(FDR = p.adjust(pValue, method = "BH"))

# Reduce results down to one representative gene set per cluster
reduced_ora_results <- full_ora_results %>%
  filter(geneSet == top_gene_set)

# Plot bar graphs of results
reduced_ora_results_split <- reduced_ora_results %>%
  group_by(ref_group, test_group) %>%
  group_split()

ora_plots <- lapply(reduced_ora_results_split, function(x) {
  x %>%
    mutate(
      log2EnrichmentRatio = log2(enrichmentRatio),
      enrichment = case_when(log2EnrichmentRatio > 0 ~ "Enriched", .default = "Depleted")
    ) %>%
    ggplot(aes(x = reorder(description, log2EnrichmentRatio), y = log2EnrichmentRatio, fill = enrichment)) +
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
    scale_fill_manual(name = "Enrichment", values = c(Enriched = "orange", Depleted = "royalblue2")) +
    coord_flip()
})
names(ora_plots) <- sapply(reduced_ora_results_split, function(x) {
  ref_group <- x$ref_group[1]
  test_group <- x$test_group[1]
  paste0(test_group, "_vs_", ref_group)
})

# Save plots to file
dir.create("plots")
for (n in names(ora_plots)) {
  p <- ora_plots[[n]]
  ggsave(paste0("plots/", cohort_id, ".ora.", n, ".png"), p)
}

# Save results to file
write_csv(full_ora_results, paste(cohort_id, "ora.full.csv", sep = "."))
write_csv(reduced_ora_results, paste(cohort_id, "ora.full.reduced.csv", sep = "."))
