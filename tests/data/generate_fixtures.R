# Run in the scrnavigator-nf/ folder
# Requires: nextflow singularity loaded; network access (copyq or equivalent).
#
# Usage (inside the annotate singularity container):
#   SIF=/scratch/er01/fj9712/.nextflow/singularity/sydneyinformaticshub-scrnavigator-nf-annotate.img
#   singularity exec $SIF Rscript --vanilla tests/data/generate_fixtures.R
#
# NOTE: Section 4 (mouse per-sample fixtures) must be generated separately in a
# Bioconductor R session with MouseGastrulationData installed — it does NOT run
# inside the annotate container. See tests/data/make_mouse_samples.R.

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
integrated_rds <- "tests/data/rds/cohort.integrated.test.rds"
if (!file.exists(integrated_rds)) {
  message(paste0(integrated_rds, " not found, creating.."))
  source("tests/data/subset_integrated.R")
} else {
  message(paste0(integrated_rds, " already exists, skipping generation."))
}

# ---------------------------------------------------------------------------
# 3. Clustered RDS subset for ANNOTATE_CUSTOM tests
# ---------------------------------------------------------------------------
clustered_rds <- "tests/data/rds/cohort.integrated.clustered.test.rds"
if (!file.exists(clustered_rds)) {
  message(paste0(clustered_rds, " not found, creating.."))
  source("tests/data/subset_clustered.R")
} else {
  message(paste0(clustered_rds, " already exists, skipping generation."))
}

# ---------------------------------------------------------------------------
# 4. Mouse per-sample RDS fixtures (3 ctrl + 3 stim) — E-MTAB-7324
#    Requires MouseGastrulationData (Bioconductor). Run separately:
#      Rscript --vanilla tests/data/make_mouse_samples.R
# ---------------------------------------------------------------------------
mouse_samples <- file.path("tests/data/rds", c(
  "mouse_ctrl1.rds", "mouse_ctrl2.rds", "mouse_ctrl3.rds",
  "mouse_stim1.rds", "mouse_stim2.rds", "mouse_stim3.rds"
))
if (!all(file.exists(mouse_samples))) {
  missing <- basename(mouse_samples[!file.exists(mouse_samples)])
  message(sprintf(
    "Mouse sample fixtures not found (%s). Generate with:\n  Rscript --vanilla tests/data/make_mouse_samples.R",
    paste(missing, collapse = ", ")
  ))
} else {
  message("Mouse sample fixtures already exist, skipping generation.")
}
