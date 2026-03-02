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
    "SingleR",
    "EnsDb.Hsapiens.v86",
    "celldex",
    "scuttle"
  ),
  site_repository = "https://cran.csiro.au/",
  ask = FALSE
)