# How to run the test suite

## Setup

Clone the repository and change into a **sibling directory**. Commands must be run from outside `projectDir` to avoid cluttering it with downloaded files and cache directories:

```bash
cd /scratch/${PROJECT}/${USER}
git clone https://github.com/Sydney-Informatics-Hub/scrnavigator-nf.git
```

### Generate the test EnsDb fixture

This step is required once before running any process-level tests that use a real EnsDb sqlite. It downloads the full human EnsDb v113 via `download_ensdb.R` inside the pipeline container.

Start an interactive job on `copyq` (has network access for pulling containers and downloading EnsDb):

```bash
qsub -I -P ${PROJECT} -q copyq -l walltime=10:00:00,mem=4GB,jobfs=200GB,storage=scratch/${PROJECT}+gdata/if89
```

Once in the session:

```bash
module load nextflow singularity
```

### Check whether the container image is already available

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

### Download EnsDb v113

Run `download_ensdb.R` inside the container. The script downloads the full human EnsDb v113 and writes `EnsDb_hsapiens_v113.sqlite` to the current directory:

```bash
singularity exec \
  --bind /scratch/${PROJECT}/${USER} \
  $SIF \
  Rscript --vanilla scrnavigator-nf/bin/download_ensdb.R 'human' 'v113' .cache
```

Move the resulting file into the test data directory:

```bash
mv EnsDb_hsapiens_v113.sqlite scrnavigator-nf/tests/data/EnsDb_hsapiens_v113.sqlite
```

This file is excluded from version control by the `*.sqlite` rule already in `.gitignore`.

To debug or inspect intermediate objects interactively, open a shell inside the container then launch R:

```bash
singularity shell \
  --bind $PWD \
  $SINGULARITY_CACHEDIR/sydneyinformaticshub-scrnavigator-nf-annotate.img
# Inside the container:
R
```

## Running tests

Tests must be run from inside `projectDir` so nf-test can find `nf-test.config`:

```bash
cd ../scrnavigator-nf
```

### Process-level tests

Require containers and the EnsDb fixture from the step above. Submit via an interactive `copyq` session or a PBS job with singularity loaded.

```bash
nf-test test tests/modules/ --profile test,singularity
```

To run a single module test:

```bash
nf-test test tests/modules/preprocess_rds.nf.test --profile test,singularity
```
