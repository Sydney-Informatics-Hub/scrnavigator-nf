#!/usr/bin/env -S Rscript --vanilla

# Install all required packages
install.packages(
  c(
    "remotes"
  ),
  repos = "https://cran.csiro.au/"
)
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")