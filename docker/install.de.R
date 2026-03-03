#!/usr/bin/env -S Rscript --vanilla

# Install all required packages
BiocManager::install(
  c(
    "DESeq2"
  ),
  site_repository = "https://cran.csiro.au/",
  ask = FALSE
)