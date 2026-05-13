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

Process-level tests run against small subset datasets committed to the repo under `tests/data/`. Run the fixture script once inside the container to generate everything needed:

```bash
singularity exec \
  $SIF \
  Rscript --vanilla scrnavigator-nf/tests/data/generate_fixtures.R
```

This creates:
- `tests/data/EnsDb_hsapiens_v113.sqlite` (~540 MB) — human Ensembl annotation database, excluded from version control
- `tests/data/rds/cohort.integrated.clustered.test.rds` — derived from the committed integrated fixture by applying Seurat clustering

Each step is skipped if the file already exists.

## Regenerating committed fixtures

The RDS fixtures below are committed to the repo and do not need to be regenerated in normal use. Follow these steps only if the pipeline's object structure has changed.

### Cell cycle RDS (`cohort.integrated.test.rds`)

Create `samplesheet.csv`:

```
sample,rds,condition,res
sample_1_ctrl,../scrnavigator-nf/tests/data/rds/sample_1.ctrl.Rds,ctrl,1
sample_3_ctrl,../scrnavigator-nf/tests/data/rds/sample_3.ctrl.Rds,ctrl,1
sample_2_stim,../scrnavigator-nf/tests/data/rds/sample_2.stim.Rds,stim,1
sample_4_stim,../scrnavigator-nf/tests/data/rds/sample_4.stim.Rds,stim,1
```

Run the pipeline through integration:

```bash
nextflow run ../scrnavigator-nf \
  --input samplesheet.csv \
  --species human \
  -resume
```

Pass the integrated RDS output path to `subset_integrated.R` — it saves the fixture to `tests/data/rds/` automatically:

```bash
singularity exec \
  $SIF \
  Rscript --vanilla scrnavigator-nf/tests/data/subset_integrated.R \
  /path/to/cohort.integrated.rds
```

Then re-run `generate_fixtures.R` to rebuild the clustered RDS from the updated integrated fixture.

### Clustered RDS (`cohort.integrated.clustered.test.rds`)

Derived from the cell cycle RDS above. Re-run `generate_fixtures.R` after updating `cohort.integrated.test.rds`.

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

## Troubleshooting

### File not found errors during nf-test

If a test fails with a `... file not found` or `checkIfExists` error, the required fixture has not been generated yet. Run the fixture script inside the container to create it:

```bash
singularity exec \
  $SIF \
  Rscript --vanilla tests/data/generate_fixtures.R
```

This generates all fixtures in order and skips any that already exist. Re-run the failing test afterwards.

## Debugging

To inspect intermediate objects interactively, open a shell inside the container then launch R:

```bash
singularity shell \
  --bind $PWD \
  $SINGULARITY_CACHEDIR/sydneyinformaticshub-scrnavigator-nf-annotate.img

# Inside the container:
R
```
