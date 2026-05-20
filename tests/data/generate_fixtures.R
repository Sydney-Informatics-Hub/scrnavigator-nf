# Run in the scrnavigator-nf/ folder
# Requires: nextflow singularity loaded; network access (copyq or equivalent).
# Also requires: tests/data/rds/cohort.integrated.test.rds
#
# Usage (inside the singularity container):
#   singularity exec ${SINGULARITY_CACHEDIR}/sydneyinformaticshub-scrnavigator-nf-annotate.img Rscript --vanilla tests/data/generate_fixtures.R

here::i_am("tests/data/generate_fixtures.R")

# ---------------------------------------------------------------------------
# 1. EnsDb v113 sqlite
# ---------------------------------------------------------------------------
ensdb_dest <- here::here("tests/data/EnsDb_hsapiens_v113.sqlite")
dl_script <- here::here("bin/download_ensdb.R")
cache_dir <- tempfile()
dir.create(cache_dir)

if (!file.exists(ensdb_dest)) {
  system2("Rscript", args = c(dl_script, "human", "v113", cache_dir))
  file.rename("EnsDb_hsapiens_v113.sqlite", ensdb_dest)
  message(sprintf("EnsDb moved to: %s (%.0f MB)", ensdb_dest, file.info(ensdb_dest)$size / 1e6))
} else {
  message("EnsDb fixture already exists, skipping download.")
}

unlink(cache_dir, recursive = TRUE)

# ---------------------------------------------------------------------------
# 2. Clustered RDS subset for ANNOTATE_CUSTOM tests
# ---------------------------------------------------------------------------
clustered_rds <- here::here("tests/data/rds/cohort.integrated.clustered.test.rds")
subset_script <- here::here("tests/data/subset_clustered.R")
if (!file.exists(clustered_rds)) {
  message(paste0(clustered_rds, " not found, creating.."))
  source(subset_script)
} else {
  message(paste0(clustered_rds, " already exists, skipping generation."))
}
