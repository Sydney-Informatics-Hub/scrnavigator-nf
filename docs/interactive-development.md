# Interactive development

## Inspecting pipeline objects in a container

To explore intermediate Seurat objects or debug R code interactively, you can run:

```bash
singularity exec \
  --bind $PWD \
  ${NXF_SINGULARITY_CACHEDIR}/sydneyinformaticshub-scrnavigator-nf-<VARIANT>.img \
  R
```

The placeholder `<VARIANT>` will be one of the following:

- `base`
- `annotate`
- `cluster`
- `doublet`
- `de`
- `fea`
- `report`

These each correspond to one of the several R containers that have been built for this pipeline and are hosted on the [`sydneyinformaticshub` DockerHub repository](https://hub.docker.com/u/sydneyinformaticshub). For example, the image `sydneyinformaticshub/scrnavigator-nf-base` will be pulled by the pipeline to `${NXF_SINGULARITY_CACHEDIR}/sydneyinformaticshub-scrnavigator-nf-base.img`.

The command above will open an R shell that you can use to load any RDS produced by the pipeline:

```r
obj <- readRDS("results/integration/cohort.integrated.rds")
```

See [testing.md](testing.md) for how to resolve and/or set `${NXF_SINGULARITY_CACHEDIR}`.
