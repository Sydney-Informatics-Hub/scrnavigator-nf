# How to run the test suite

## Setup

Clone the repository and change into a **sibling directory**. Commands must be run from outside `projectDir` to avoid cluttering it with downloaded files and cache directories:

```bash
cd /scratch/${PROJECT}/${USER}
git clone https://github.com/Sydney-Informatics-Hub/scrnavigator-nf.git
mkdir run && cd run   # work from here, never from inside scrnavigator-nf/
```

All paths in the commands below assume this layout:

```
/scratch/${PROJECT}/${USER}/
├── run/                  ← run commands from here
└── scrnavigator-nf/      ← projectDir, do not cd into
```

## Generate the test EnsDb fixture

This step is required once before running any process-level tests that use a real EnsDb sqlite. It downloads the full human EnsDb and subsets it to only the genes present in the test data.

Start an interactive job on `copyq` (has network access for pulling containers and downloading EnsDb):

```bash
qsub -I -P ${PROJECT} -q copyq -l walltime=10:00:00,mem=4GB,storage=scratch/${PROJECT}+gdata/if89
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

**If the image is listed**, use it directly:

```bash
SIF="${SINGULARITY_CACHEDIR}/sydneyinformaticshub-scrnavigator-nf-annotate.img"
singularity exec --bind $PWD $SIF Rscript scrnavigator-nf/tests/data/subset_sqlite.R
```

**If the image is not listed**, pull it first then run:

```bash
singularity pull \
  --dir ${SINGULARITY_CACHEDIR} \
  docker://sydneyinformaticshub/scrnavigator-nf-annotate

SIF="${SINGULARITY_CACHEDIR}/sydneyinformaticshub-scrnavigator-nf-annotate.img"
singularity exec --bind $PWD $SIF Rscript scrnavigator-nf/tests/data/subset_sqlite.R
```

To run the script interactively (e.g. to debug or inspect intermediate objects), open a shell inside the container then launch R:

```bash
singularity shell --bind $PWD $SIF
# inside the container open an interactive terminal:
R
```

This produces `tests/data/EnsDb_hsapiens_test.sqlite`. This file is excluded from version control by the `*.sqlite` rule already in `.gitignore`.
