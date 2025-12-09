# Get commandline arguments
args <- commandArgs(trailingOnly = TRUE)

sample_id <- args[1]
rds_path <- args[2]
is_ensembl <- as.logical(args[3])
if (is.na(is_ensembl)) {
    stop("Error: Invalid argument for 'is_ensembl'")
}

# TODO