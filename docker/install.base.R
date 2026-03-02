#!/usr/bin/env -S Rscript --vanilla

# Install all required packages
install.packages(
  c(
    "tidyverse",
    "Seurat",
    "ggplot2",
    "ggrepel",
    "hdf5r"
  ),
  repos = "https://cran.csiro.au/"
)