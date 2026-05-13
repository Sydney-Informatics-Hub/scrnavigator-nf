# How to run the test suite

## Setup

Clone the repository:

```bash
cd /scratch/${PROJECT}/${USER}
git clone https://github.com/Sydney-Informatics-Hub/scrnavigator-nf.git
```

### Start an interactive session

Process-level tests need a Singularity container and (for fixture generation) network access. Submit an interactive job on `copyq`:

```bash
qsub -I -P ${PROJECT} -q copyq -l walltime=10:00:00,mem=4GB,jobfs=200GB,storage=scratch/${PROJECT}+gdata/if89
```

Once in the session:

```bash
module load nextflow singularity
```

### Resolve the container image

The pipeline caches Singularity images under your Nextflow singularity cache directory. Check whether the annotate image has been pulled already:

```bash
ls $SINGULARITY_CACHEDIR | grep scrnavigator-nf-annotate
```

**If the image is listed**, set `SIF` to its path:

```bash
SIF=${SINGULARITY_CACHEDIR}/sydneyinformaticshub-scrnavigator-nf-annotate.img
```

**If the image is not listed**, pull it first:

```bash
singularity pull \
  --dir ${SINGULARITY_CACHEDIR} \
  docker://sydneyinformaticshub/scrnavigator-nf-annotate

SIF=$(ls $SINGULARITY_CACHEDIR | grep scrnavigator-nf-annotate | head -1)
SIF="${SINGULARITY_CACHEDIR}/${SIF}"
```

## Generate the test fixtures

Process-level tests run against small subset datasets that are committed to the repo under `tests/data/`. Subsets keep the test suite fast and let it run on modest hardware while still exercising real container behaviour, R code paths, and file I/O.

The exception is the Ensembl annotation database (`EnsDb_hsapiens_v113.sqlite`, ~1 GB), which is too large for version control and must be generated locally before running tests that depend on it.

### EnsDb v113 sqlite

Run `download_ensdb.R` inside the container. It downloads the full human EnsDb v113 and writes `EnsDb_hsapiens_v113.sqlite` to the current directory:

```bash
singularity exec \
  $SIF \
  Rscript --vanilla scrnavigator-nf/bin/download_ensdb.R 'human' 'v113' .cache
```

Move the resulting file into the test data directory:

```bash
mv EnsDb_hsapiens_v113.sqlite scrnavigator-nf/tests/data/EnsDb_hsapiens_v113.sqlite
```

This file is excluded from version control by the `*.sqlite` rule in `.gitignore`.

### Cell cycle RDS (already committed)

`tests/data/rds/cohort.integrated.test.rds` is committed to the repo and ready to use. The instructions below are only needed if you want to regenerate it (e.g. after pipeline changes that affect the integrated object structure).

Create `samplesheet.csv`:

```
sample,rds,condition,res
sample_1_ctrl,../scrnavigator-nf/tests/data/rds/sample_1.ctrl.Rds,ctrl,1
sample_3_ctrl,../scrnavigator-nf/tests/data/rds/sample_3.ctrl.Rds,ctrl,1
sample_2_stim,../scrnavigator-nf/tests/data/rds/sample_2.stim.Rds,stim,1
sample_4_stim,../scrnavigator-nf/tests/data/rds/sample_4.stim.Rds,stim,1
```

Then run the pipeline through integration:

```bash
nextflow run ../scrnavigator-nf \
  --input samplesheet.csv \
  --species human \
  -resume
```

Copy the integrated RDS output to `tests/data/rds/cohort.integrated.rds`, then subset it to a small test fixture (~10 MB):

```bash
singularity exec \
  $SIF \
  Rscript --vanilla scrnavigator-nf/tests/data/subset_integrated.R
```

This produces `tests/data/rds/cohort.integrated.test.rds`.

## Running tests

Tests must be run from inside `projectDir` so nf-test can find `nf-test.config`:

```bash
cd scrnavigator-nf
```

### Process-level tests

Require containers and the EnsDb fixture from the step above. Submit via an interactive `copyq` session (recommended) or a PBS job with singularity loaded.

```bash
nf-test test tests/modules/ --profile test,singularity
```

To run a single module test:

```bash
nf-test test tests/modules/preprocess_rds.nf.test --profile test,singularity
```

## Debugging

To inspect intermediate objects interactively, open a shell inside the container then launch R:

```bash
singularity shell \
  --bind $PWD \
  $SINGULARITY_CACHEDIR/sydneyinformaticshub-scrnavigator-nf-annotate.img

# Inside the container:
R
```
