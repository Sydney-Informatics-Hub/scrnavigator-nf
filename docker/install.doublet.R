#!/usr/bin/env -S Rscript --vanilla

# Install all required packages
install.packages(
  c(
    "clustree"
  ),
  repos = "https://cran.csiro.au/"
)
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")