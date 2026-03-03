#!/usr/bin/env -S Rscript --vanilla

# Install all required packages
install.packages(
  c(
    "tidyverse",
    "Seurat",
    "ggplot2",
    "ggrepel",
    "hdf5r",
    "BiocManager",
    "remotes"
  ),
  repos = "https://cran.csiro.au/"
)
remotes::install_github("immunogenomics/presto")
remotes::install_github("bnprks/BPCells/r")
BiocManager::install(
  c(
    "glmGamPoi"
  ),
  site_repository = "https://cran.csiro.au/",
  ask = FALSE
)