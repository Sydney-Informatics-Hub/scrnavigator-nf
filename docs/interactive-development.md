# Interactive development

## Inspecting pipeline objects in a container

To explore intermediate Seurat objects or debug R code interactively, open a shell inside the annotate container and launch R from there:

```bash
singularity shell \
  --bind $PWD \
  $SINGULARITY_CACHEDIR/sydneyinformaticshub-scrnavigator-nf-annotate.img

# Inside the container:
R
```

From the R session you can load any RDS produced by the pipeline:

```r
obj <- readRDS("results/integration/cohort.integrated.rds")
```

See [testing.md](testing.md) for how to resolve `$SINGULARITY_CACHEDIR` and set `SIF` if you haven't done so already.
