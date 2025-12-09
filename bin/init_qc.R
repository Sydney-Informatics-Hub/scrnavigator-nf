library(tidyverse)

# Get commandline arguments
args <- commandArgs(trailingOnly = TRUE)

sample_id <- args[1]
rds_path <- args[2]
params_file <- args[3]
template_file <- args[4]

# Read in parameters
params <- read_csv(params_file)

# TODO

rmarkdown::render(template_file, params = list(qc_files = "qc_files"))