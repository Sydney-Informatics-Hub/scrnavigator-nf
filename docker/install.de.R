#!/usr/bin/env -S Rscript --vanilla

# Install all required packages
install.packages(
  c(
    "BiocManager"
  ),
  repos = "https://cran.csiro.au/"
)
BiocManager::install(
  c(
    "DESeq2"
  ),
  site_repository = "https://cran.csiro.au/",
  ask = FALSE
)