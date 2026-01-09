#!/usr/bin/env -S Rscript --vanilla
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggrepel)

# Get commandline arguments
args <- commandArgs(trailingOnly = TRUE)

cohort_id <- args[1]
gene_symbols_file <- args[2]
params_file <- args[3]
rds_paths <- args[4:length(args)]

# Get parameters
params <- read_csv(params_file)
p_val_cutoff <- as.numeric(params$value[match("p_val_cutoff", params$param)])
fc_cutoff <- as.numeric(params$value[match("fc_cutoff", params$param)])

# Create output directory
dir.create("results")

# Get gene symbols dataframe
gene_symbols <- readRDS(gene_symbols_file)

# Read in DE results and update significance values
all_de_results <- lapply(rds_paths, readRDS) %>%
  bind_rows() %>%
  mutate(
    p_val_adj_all = p.adjust(p_val, method = "bonferroni"),
    neg_log10_pval_adj_all = -log10(p_val_adj_all),
    sig_p = p_val_adj_all < p_val_cutoff
  )

if (!is.na(fc_cutoff)) {
  all_de_results <- all_de_results %>%
    mutate(
      sig = sig_p & (
        avg_log2FC <= -log2(fc_cutoff) |
          avg_log2FC >= log2(fc_cutoff)
      )
    )
} else {
  all_de_results <- all_de_results %>%
    mutate(sig = sig_p)
}

all_de_results <- all_de_results %>%
  mutate(
    sig = case_when(
      sig ~ "sig",
      .default = "ns"
    ),
    sig_p = case_when(
      sig_p ~ "sig",
      .default = "ns"
    )
  )

# Add gene symbols if using Ensembl IDs
all_de_results$Gene <- rownames(all_de_results)
using_ens_ids <- all(startsWith(all_de_results$Gene, "ENS"))
if (using_ens_ids) {
  idx <- match(all_de_results$Gene, gene_symbols$gene_versions)
  all_de_results$Gene.Symbol <- gene_symbols$gene_symbols[idx]
} else {
  all_de_results$Gene.Symbol <- all_de_results$Gene
}

# Plot p-value and log-fold-change distributions
all_de_results <- all_de_results %>%
  unite(comparison, test_group, ref_group, sep = " vs. ", remove = FALSE)

p_p_val_distributions <- all_de_results %>%
  drop_na(p_val) %>%
  ggplot(aes(x = p_val)) +
  geom_histogram(binwidth = 0.01) +
  facet_wrap(facets = ~ comparison, ncol = 2) +
  theme_light() +
  annotation_logticks(sides = "l", colour = "lightgrey")
ggsave(paste0("results/", cohort_id, ".p_val_dist.png"), p_p_val_distributions)

p_log_fc_distributions <- all_de_results %>%
  drop_na(avg_log2FC, sig_p) %>%
  ggplot(aes(x = avg_log2FC, fill = sig_p)) +
  geom_histogram(binwidth = 0.2, position = "identity", alpha = 0.5) +
  geom_vline(xintercept = c(-log2(fc_cutoff), log2(fc_cutoff))) +
  scale_y_log10() +
  facet_wrap(facets = ~ comparison, ncol = 1) +
  coord_cartesian(xlim = c(-5, 5)) +
  theme_light() +
  annotation_logticks(sides = "l", colour = "lightgrey")
ggsave(paste0("results/", cohort_id, ".log_fc_dist.png"), p_log_fc_distributions)

# Plot volcano plots
top10 <- all_de_results %>%
  group_by(comparison) %>%
  arrange(p_val_adj_all) %>%
  dplyr::filter(sig == "sig") %>%
  slice_head(n = 10) %>%
  ungroup

p_volcano_plots <- all_de_results %>%
  drop_na(avg_log2FC, neg_log10_pval_adj_all, sig, Gene.Symbol) %>%
  ggplot(aes(x = avg_log2FC, y = neg_log10_pval_adj_all, colour = sig)) +
  geom_point() +
  geom_hline(yintercept = -log10(p_val_cutoff)) +
  scale_colour_manual(breaks = c("ns", "sig"), values = c("#bababa", "#ca0020")) +
  geom_text_repel(data = top10, aes(label = Gene.Symbol)) +
  facet_wrap(facets = ~ comparison, ncol = 1) +
  xlab(expression(log[2] * FC)) + ylab(expression(-log[10] * p.value.bonf)) +
  theme_light()
if (!is.na(fc_cutoff)) {
  p_volcano_plots <- p_volcano_plots + geom_vline(xintercept = c(-log2(fc_cutoff), log2(fc_cutoff)))
}
ggsave(paste0("results/", cohort_id, ".volcano.png"), p_volcano_plots)

# Get all significantly DE genes
all_de_results %>%
  dplyr::filter(sig == "sig") %>%
  write_csv(paste0("results/", cohort_id, ".significant_de_genes.all_tests.csv"))

# Save differential expression data to file
write_csv(all_de_results, paste(cohort_id, "de.full.csv", sep = "."))
saveRDS(all_de_results, paste(cohort_id, "de.full.Rds", sep = "."))
