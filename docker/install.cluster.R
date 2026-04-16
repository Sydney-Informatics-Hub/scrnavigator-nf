#!/usr/bin/env -S Rscript --vanilla

# Install all required packages
install.packages(
  c(
    "clustree",
    "leidenbase"
  ),
  repos = "https://cran.csiro.au/"
)