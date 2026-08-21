#!/usr/bin/env -S Rscript --vanilla

# Install all required packages
BiocManager::install(
  c(
    "SingleR",
    "EnsDb.Hsapiens.v86",
    "EnsDb.Mmusculus.v79",
    "celldex",
    "scuttle"
  ),
  site_repository = "https://cran.csiro.au/",
  ask = FALSE
)