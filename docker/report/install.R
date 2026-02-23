#!/usr/bin/env -S Rscript --vanilla

# Install all required packages
install.packages(
  c(
    "tidyverse",
    "DT"
  ),
  repos = "https://cran.csiro.au/"
)
