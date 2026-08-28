# Run in the scrnavigator-nf/ folder
# Requires: nextflow singularity loaded; network access (copyq or equivalent).
# Also requires: tests/data/rds/cohort.integrated.test.rds
#
# Usage (inside the annotate singularity container):
#   SIF=/scratch/er01/fj9712/.nextflow/singularity/sydneyinformaticshub-scrnavigator-nf-annotate.img
#   singularity exec $SIF Rscript --vanilla tests/data/generate_fixtures.R
#
# NOTE: Section 4 (mouse per-sample fixtures) must be generated separately in a
# Bioconductor R session with MouseGastrulationData installed — it does NOT run
# inside the annotate container. See tests/data/make_mouse_samples.R.

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
