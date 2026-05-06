# Run in the scrnavigator-nf/ folder
# Requires: nextflow singularity loaded; network access (copyq or equivalent).
#
# Usage (inside the singularity container):
#   singularity exec ${SINGULARITY_CACHEDIR}/sydneyinformaticshub-scrnavigator-nf-annotate.img Rscript --vanilla tests/data/generate_fixtures.R

# ---------------------------------------------------------------------------
# 1. EnsDb v113 sqlite
# ---------------------------------------------------------------------------
ensdb_dest <- "tests/data/EnsDb_hsapiens_v113.sqlite"

if (!file.exists(ensdb_dest)) {
  system2("Rscript", args = c("bin/download_ensdb.R", "human", "v113", ".cache"))
  file.rename("EnsDb_hsapiens_v113.sqlite", ensdb_dest)
  message(sprintf("EnsDb moved to: %s (%.0f MB)", ensdb_dest, file.info(ensdb_dest)$size / 1e6))
} else {
  message("EnsDb fixture already exists, skipping download.")
}

# ---------------------------------------------------------------------------
# 2. Integrated RDS subset for ANNOTATE_CELL_CYCLE tests
# ---------------------------------------------------------------------------
source("tests/data/subset_integrated.R")