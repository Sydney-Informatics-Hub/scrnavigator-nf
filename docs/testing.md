# Testing

The following instructions demonstrate how to set up and run the pipeline tests via nf-test. Note that some of the instructions are specific to running these tests on the Australian national HPC "Gadi", hosted by the National Computing Infrastructure (NCI); these steps have been marked as such. These instructions also assume you are using Singularity as your container software.

## Setup

Clone the repository:

```bash
cd /scratch/${PROJECT}/${USER}
git clone https://github.com/Sydney-Informatics-Hub/scrnavigator-nf.git
```

### Gadi: start an interactive session

Process-level tests need a Singularity container and (for fixture generation) network access. For running these tests on Gadi, you will need to submit an interactive job on `copyq`:

```bash
qsub -I -P ${PROJECT} -q copyq -l walltime=10:00:00,mem=8GB,jobfs=200GB,storage=scratch/${PROJECT}
```

Once in the session:

```bash
module load nextflow singularity
```

### Resolve the container image

The pipeline caches Singularity images under your Nextflow singularity cache directory. By default, this is located at `<NXF_LAUNCH_DIR>/work/singularity`, where `<NXF_LAUNCH_DIR>` is the launch directory where Nextflow is run from. It is recommended that you set this to a different directory of your choosing, especially if you are running Nextflow on a shared system like an HPC. You can set the Nextflow singularity cache directory by exporting the environment variable `NXF_SINGULARITY_CACHEDIR`.

Generating the test fixures requires the `sydneyinformaticshub/scrnavigator-nf-annotate` image. Check whether the annotate image has been pulled already (e.g. if you have previously run the pipeline):

```bash
SIF=${NXF_SINGULARITY_CACHEDIR}/sydneyinformaticshub-scrnavigator-nf-annotate.img
ls ${SIF}
```

If the image is **not listed**, pull it now:

```bash
singularity pull \
  ${SIF} \
  docker://sydneyinformaticshub/scrnavigator-nf-annotate
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

It will also re-create the following committed fixture if it is missing:
- `tests/data/rds/cohort.integrated.clustered.test.rds` — derived from the committed integrated fixture by applying Seurat clustering

Each step is skipped if the file already exists.

## Regenerating committed fixtures

The following RDS fixtures are committed to the repo and do not need to be regenerated in normal use:

- `tests/data/cohort.integrated.test.rds`
- `tests/data/cohort.integrated.clustered.test.rds`

Follow these steps only if the pipeline's object structure has changed and new versions of these files are required.

These steps all assume you are running from the project directory's **parent directory**.

First, delete the existing fixtures:

```bash
rm scrnavigator-nf/tests/data/cohort.integrated.test.rds
rm scrnavigator-nf/tests/data/cohort.integrated.clustered.test.rds
```

Run the pipeline through integration using the samplesheet provided in `tests/data/samplesheet.generate_fixtures.csv`:

```bash
# Assuming running from the project directory's parent
nextflow run scrnavigator-nf \
  --input scrnavigator-nf/tests/data/samplesheet.generate_fixtures.csv \
  --species human \
  --no_analysis true \
  -resume
```

Pass the integrated RDS output path to `subset_integrated.R` — it saves the fixture to `tests/data/rds/` automatically:

```bash
singularity exec \
  $SIF \
  Rscript --vanilla scrnavigator-nf/tests/data/subset_integrated.R \
  scrnavigator-nf/tests/data/cohort.integrated.test.rds
```

Then re-run `generate_fixtures.R` [as described above](#generate-the-test-fixtures) to rebuild the clustered RDS from the updated integrated fixture.

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
