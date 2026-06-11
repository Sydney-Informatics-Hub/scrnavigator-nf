# scrnavigator-nf

## Description

This Nextflow pipeline provides a standardised method for processing and analysing single cell RNA sequencing (scRNAseq) data. It builds upon the [scRNAvigator notebooks](https://github.com/Sydney-Informatics-Hub/scrna-analysis) - a series of interactive Quarto notebooks that step the user through quality control, dataset integration, cell annotation, pseudobulking, differential expression and functional enrichment analysis (FEA). In this Nextflow pipeline, these same steps have been organised into a semi-automated workflow that can be easily scaled up and take advantage of parallelisation and high performance computing infrastructures. The workflow is designed to be run in an iterative way, starting with quality control and filtering, then integration and cell annotation, and finally pseudobulking, differential expression and FEA. This allows the user to inspect vital quality metrics at each step and select filtering thresholds, annotations and sample groupings required for the downstream stages.

## User guide

### Pipeline setup

Download the repo to your launch directory.

```bash
git clone git@github.com:Sydney-Informatics-Hub/scrnavigator-nf.git
```

### Input data

The initial input data that this pipeline requires is one or more RDS files, each containing a Seurat object with the count matrix data for a single sample. [Seurat is an R package](https://satijalab.org/seurat/) that provides a framework for handling and processing single cell RNAseq data.

The initial data is expected to have already been pre-processed to remove background noise, i.e. ambient RNA signals. This is done automatically by 10X's `cellranger` software for their platform, and similar filtering can be achieved with tools like [`cellbender`](https://github.com/broadinstitute/CellBender) for other datasets.

The Seurat data object should at the very minimum contain the count matrix data for your sample in its `RNA` assay, along with either Ensembl gene IDs or HGNC gene symbols as the `RNA` assay's row names. The Seurat object's `meta.data` field should also be present with a row for every cell in the sample.

This pipeline was designed to work with the filtered output from the [nf-core `scrnaseq` Nextflow pipeline](https://nf-co.re/scrnaseq/). The `scrnaseq` pipeline can handle various scRNAseq datasets, including those from the 10X platform, and will generate RDS files containing the Seurat data required for this pipeline. We highly recommend using the `scrnaseq` pipeline for the initial alignment, counting and pre-processing of your data prior to using this workflow.

### Example samplesheet

## Parameters

See [docs/params.md](docs/params.md) for the full parameter reference.

## Usage

> **Note:** The pipeline is currently under active development and usage may change rapidly.

Run each step sequentially, inspecting the interactive report (`results/report/index.html`) before proceeding to the next. Outputs from each step are cached, so re-running with `-resume` skips completed work.

### Quality control and filtering

Start with an unfiltered QC pass to understand the quality characteristics of each sample before applying any cell filters. The `--qc_only` flag runs only the QC subworkflow, producing per-sample metrics to guide threshold selection.

**1. Create a samplesheet** with at minimum a sample name and path to a preprocessed Seurat RDS file. Additional metadata columns can be added, in this example, the donor `sex`:

```console
sample,rds,sex
Donor_1,/g/data/er01/test-data/single-cell/human/rds/Donor1_filtered_matrix.seurat.rds,male
Donor_2,/g/data/er01/test-data/single-cell/human/rds/Donor2_filtered_matrix.seurat.rds,male
Donor_3,/g/data/er01/test-data/single-cell/human/rds/Donor3_filtered_matrix.seurat.rds,female
Donor_4,/g/data/er01/test-data/single-cell/human/rds/Donor4_filtered_matrix.seurat.rds,female
```

**2. Run the initial QC pass:**

```bash
nextflow run scrnavigator-nf \
    --input path/to/samplesheet.csv \
    --outdir results \
    --species human \
    --qc_only
```

Outputs:

```
results/
├── qc
│   ├── Donor_1
│   │   ├── cluster    # Clustering plots and barcode cell assignments across resolutions
│   │   ├── filter     # Low-quality cell filtering metrics
│   │   └── preprocess
│   ├── Donor_2
│   │   └── ...
├── report
│   └── report         # Interactive HTML report to explore QC and clustering results
└── run_info
```

**3. Inspect the report.** Open `results/report/index.html` and identify per-sample filtering thresholds by examining the distributions of:

- **nCount_RNA** - total RNA counts per cell
- **nFeature_RNA** - unique genes per cell
- **percent.mt** - mitochondrial gene fraction

**4. Update your samplesheet** with per-sample filtering thresholds. Leave any bound blank if no threshold applies for that sample:

```console
sample,rds,sex,min_ncount,max_ncount,min_nfeature,max_nfeature,min_mt_pct,max_mt_pct
Donor_1,/g/data/er01/test-data/single-cell/human/rds/Donor1_filtered_matrix.seurat.rds,male,3000,12000,2000,20000,,10
Donor_2,/g/data/er01/test-data/single-cell/human/rds/Donor2_filtered_matrix.seurat.rds,male,3000,12000,2000,20000,,10
Donor_3,/g/data/er01/test-data/single-cell/human/rds/Donor3_filtered_matrix.seurat.rds,female,3000,12000,2000,20000,,10
```

**5. Re-run with thresholds applied,** using `-resume` to skip the already-completed unfiltered pass:

```bash
nextflow run scrnavigator-nf \
    --input path/to/samplesheet.csv \
    --outdir results \
    --species human \
    --qc_only \
    -resume
```

Once satisfied with the per-sample filtering, proceed to doublet detection.

### Doublet detection and dataset integration

Next, we proceed to doublet detection and integration. 

Doublets are droplets that captured two cells and are flagged using DoubletFinder. Once removed, samples are merged, batch-corrected using CCA integration, and clustered. This is run by dropping the `--qc-only` flag. The `--no_analysis` flag is set to prevent moving on to cell annotation and analysis before inspecting the doublet and integration QC results (command specified below).

Before running, **add a clustering resolution column (`res`) to the samplesheet** for each sample. The resolution controls the granularity of per-sample clusters used to inform doublet detection. Review the clustree plots from the QC step to choose an appropriate value per sample:

```console
sample,rds,sex,min_ncount,max_ncount,min_nfeature,max_nfeature,min_mt_pct,max_mt_pct,res
Donor_1,/g/data/er01/test-data/single-cell/human/rds/Donor1_filtered_matrix.seurat.rds,male,3000,12000,2000,20000,,10,1.2
Donor_2,/g/data/er01/test-data/single-cell/human/rds/Donor2_filtered_matrix.seurat.rds,male,3000,12000,2000,20000,,10,1
Donor_3,/g/data/er01/test-data/single-cell/human/rds/Donor3_filtered_matrix.seurat.rds,female,3000,12000,2000,20000,,10,1.2
Donor_4,/g/data/er01/test-data/single-cell/human/rds/Donor4_filtered_matrix.seurat.rds,female,3000,12000,2000,20000,,10,0.6
```

Replace `--qc_only` with `--no_analysis` and re-run. QC steps will be retrieved from cache:

```bash
nextflow run scrnavigator-nf \
    --input path/to/samplesheet.csv \
    --outdir results \
    --species human \
    --no_analysis \
    -resume
```

You should see QC steps served from cache, followed by doublet detection and integration running for the first time:

```
executor >  pbspro (4)
[39/1b16c5] QUALITY_CONTROL:PREPROCESS_RDS (3) [100%] 4 of 4, cached: 4
[81/850960] QUALITY_CONTROL:FILTER (4)         [100%] 4 of 4, cached: 4
[1b/333776] QUALITY_CONTROL:SCTRANSFORM (4)    [100%] 4 of 4, cached: 4
[fa/fd3d60] DETECT_DOUBLETS (4)                [  0%] 0 of 4
[-        ] INTEGRATE:INTEGRATION              -
[-        ] INTEGRATE:CLUSTER_INTEGRATED       -
[-        ] REPORT                             -
```

Outputs:

```
results/
├── integration                               # NEW: cohort-level integration outputs
│   ├── clustering
│   │   ├── cohort.integrated.clustered.rds  # Integrated Seurat object with cluster assignments at multiple resolutions
│   │   └── qc_results                       # Clustree, UMAP, and per-cluster feature/count plots
│   ├── cohort.integrated.rds                # Integrated Seurat object (CCA-integrated, UMAP computed, pre-clustering)
│   ├── gene_symbols.Rds                     # Dataframe mapping gene symbols to Ensembl IDs for the cohort
│   └── qc_results
│       ├── cohort.umap.integrated.png       # UMAP coloured by sample after CCA integration (batch correction applied)
│       └── cohort.umap.merged.png           # UMAP coloured by sample before integration (merged only)
├── qc
│   ├── Donor_1
│   │   ├── cluster
│   │   ├── doublets                         # NEW: doublet detection results and plots per sample
│   │   │   ├── Donor_1.doublets_detected.rds              # Putative doublets flagged, prior to removal
│   │   │   ├── Donor_1.doublets_removed.sct_clustered.rds # Seurat object with doublets removed, re-normalised and clustered
│   │   │   └── qc_results                                 # Doublet summary table and diagnostic plots for the report
│   │   │       ├── Donor_1.doublet_summary.csv            # Per-cluster doublet counts and proportions
│   │   │       ├── Donor_1.doublet_umap.1.png             # UMAP coloured by doublet status
│   │   │       ├── Donor_1.doublets_per_cluster.1.png     # Doublet counts per cluster
│   │   │       ├── Donor_1.doublet_proportions_per_cluster.1.png
│   │   │       ├── Donor_1.clustree.png                   # Resolution tree across clustering parameters
│   │   │       └── Donor_1.feature_count_plot.cluster.*.png # nCount vs nFeature per cluster, one plot per resolution
│   │   ├── filter
│   │   └── preprocess
│   ├── Donor_2
│   │   └── ...
├── report
│   └── report                               # UPDATED: Interactive HTML report
└── run_info
```

Inspect the integration UMAPs in the report to confirm that samples mix well after batch correction. If samples remain separated post-integration, consider revisiting QC thresholds or integration parameters.

### Cell annotation

Cell annotation assigns biological identity to each cell using three complementary approaches: database-driven annotation (SingleR), cell cycle phase scoring, and custom marker gene-based annotation. Cell cycle scoring is performed automatically for human samples; databased-driven annotation is also performed automatically for both human and mouse samples, using the [Human Primary Cell Atlas (HPCA)](https://www.bioconductor.org/packages/release/data/experiment/vignettes/celldex/inst/doc/userguide.html#1_Human_Primary_Cell_Atlas) and [MouseRNAseqData](https://www.bioconductor.org/packages/release/data/experiment/vignettes/celldex/inst/doc/userguide.html#2_Mouse_RNA-seq_Data) databases, respectively, both sourced from the `celldex` package. See [How-to: use a custom annotation database](#how-to-use-a-custom-annotation-database) for using an alternative reference database with SingleR.

Along with cell-level annotation, cell clusters are labelled with a cell type if that type represents more than 67% of cells in the cluster; otherwise it is marked "Ambiguous". The defaul 67% threshold can be overridden with the `--cell_type_proportion_threshold` parameter.

Drop `--no_analysis` to run annotation (and pseudobulking, described below) for the first time:

```bash
nextflow run scrnavigator-nf \
    --input path/to/samplesheet.csv \
    --outdir results \
    --species human \
    -resume
```

Outputs:

```
results/
├── annotation                               # NEW
│   ├── cell_cycle
│   │   ├── available_annotations.txt        # Lists annotation fields added (e.g. 'Phase')
│   │   ├── cohort.annotated.cell_cycle.rds  # Seurat object with Phase (G1/S/G2M) added to metadata
│   │   └── qc_results                       # UMAP coloured by cell cycle phase
│   ├── clusters
│   │   ├── cohort.annotated.clusters.rds           # Seurat object with cluster_annotation column added to metadata
│   │   ├── cohort.cluster_cell_type_assignments.csv # Mapping of cluster IDs to consensus cell type labels
│   │   └── qc_results                              # Cell type proportion dot plots and tables per cluster
│   └── db
│       ├── available_annotations.txt               # Lists SingleR annotation fields added (e.g. 'SingleR.hpca_main', 'SingleR.hpca_fine')
│       ├── cohort.annotated.database.rds            # Seurat object with SingleR cell type annotations added to metadata
│       └── qc_results                              # UMAP and cell type proportion plots for SingleR annotations
└── report
    ├── report                               # UPDATED: Interactive HTML report
    └── report.tar
```

Inspect the annotation UMAPs and cell type proportion plots in the report. Verify that cluster assignments look biologically plausible before defining comparison groups for downstream analysis.

### How-to: use a custom annotation database

By default the pipeline annotates human samples against the [Human Primary Cell Atlas](https://www.bioconductor.org/packages/release/data/experiment/vignettes/celldex/inst/doc/userguide.html#1_Human_Primary_Cell_Atlas) and mouse samples against [MouseRNAseqData](https://www.bioconductor.org/packages/release/data/experiment/vignettes/celldex/inst/doc/userguide.html#2_Mouse_RNA-seq_Data), both sourced from the `celldex` package. You can supply any alternative reference by providing a path to an RDS file via the `--annotation_db` parameter.

**Required format**

The database must be a [`SummarizedExperiment`](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html) object with:

- A `logcounts` assay containing log-normalised expression values, with genes as rows and reference cells/samples as columns.
- A `label` column in `colData` containing the cell type label for each reference cell/sample.

Row names must use the same gene identifier type as your data (Ensembl IDs or gene symbols).

**Building a custom database**

The example below constructs a minimal compliant object from a raw counts matrix and a vector of cell type labels:

```r
library(SummarizedExperiment)
library(scuttle)

# counts_matrix: genes x cells, raw counts
# cell_labels:   character vector, one label per cell

sce <- SummarizedExperiment(
  assays = list(counts = counts_matrix),
  colData = DataFrame(label = cell_labels)
)
sce <- logNormCounts(sce)   # adds the required 'logcounts' assay

saveRDS(sce, "my_annotation_db.rds")
```

Any `celldex` reference is already in this format and can be saved and reused directly:

```r
db <- celldex::HumanPrimaryCellAtlasData()
saveRDS(db, "hpca.rds")
```

**Running the pipeline with a custom database**

Pass the RDS path with `--annotation_db`. The `--species` parameter is still required for other steps in the pipeline but does not affect which database is used when `--annotation_db` is provided:

```bash
nextflow run scrnavigator-nf \
    --input path/to/samplesheet.csv \
    --outdir results \
    --species human \
    --annotation_db path/to/my_annotation_db.rds \
    -resume
```

The resulting per-cell annotations will be stored in the `SingleR.annotation` metadata column of the output Seurat object.

### Pseudobulking

Pseudobulking aggregates single-cell counts per sample per group, creating sample-level expression profiles suitable for bulk-style DE testing. Groups are defined by combinations of metadata fields (e.g., sex and cell type), and each unique combination becomes one pseudobulk sample. A minimum of 3 samples per group is recommended for reliable DE results.

**1. Inspect the annotation report.** Identify which annotation field best captures cell identity (e.g., `cluster_annotation` from the clusters step).

**2. Create a `comparisons.csv`** defining the pairwise contrasts to test:

```console
ref,test
female_B-cell,male_B-cell
female_T-cells,male_T-cells
female_NK-cell,male_NK-cell
female_Monocyte,male_Monocyte
```

**3. Re-run** with `--pseudo_groups` and `--comparisons` added:

```bash
nextflow run scrnavigator-nf \
    --input path/to/samplesheet.csv \
    --outdir results \
    --species human \
    --pseudo_groups sex,cluster_annotation \
    --comparisons path/to/comparisons.csv \
    -resume
```

Outputs:

```
results/
├── analysis
│   └── pseudobulk
│       ├── cohort.comparison_groups.txt     # Comparison groups and sample counts per group
│       └── cohort.pseudobulk.rds            # Seurat object with summed counts aggregated per sample per group
├── annotation
└── report
    ├── report                               # UPDATED: Interactive HTML report
    └── report.tar
```

Check `cohort.comparison_groups.txt` to confirm each group contains enough samples (≥ 3 recommended) before proceeding to DE analysis.

### Differential expression analysis

Pairwise DE testing is performed between conditions defined in `comparisons.csv` using DESeq2, operating on the pseudobulked samples. Results are combined across all comparisons and a Bonferroni correction is applied across all tests to control for multiple comparisons.

Outputs:

```
results/
├── analysis
│   ├── differential_expression              # NEW
│   │   ├── cohort.de.full.csv               # All comparisons combined; includes Bonferroni-adjusted p-value and significance flags
│   │   ├── cohort.de.full.Rds               # Same as above (R format)
│   │   ├── cohort.de.male_B-cell.female_B-cell.csv  # Per-comparison DE results (log2FC, p-value, pct.1/pct.2 per gene)
│   │   ├── cohort.de.male_B-cell.female_B-cell.Rds  # Same as above (R format)
│   │   ├── cohort.de.male_Monocyte.female_Monocyte.csv
│   │   ├── cohort.de.male_Monocyte.female_Monocyte.Rds
│   │   ├── cohort.de.male_NK-cell.female_NK-cell.csv
│   │   ├── cohort.de.male_NK-cell.female_NK-cell.Rds
│   │   ├── cohort.de.male_T-cells.female_T-cells.csv
│   │   ├── cohort.de.male_T-cells.female_T-cells.Rds
│   │   └── results
│   │       ├── cohort.p_val_dist.png        # P-value distribution histogram, faceted by comparison
│   │       ├── cohort.log_fc_dist.png       # Log2FC distribution coloured by significance, faceted by comparison
│   │       ├── cohort.volcano.png           # Volcano plot with top 10 significant genes labelled per comparison
│   │       └── cohort.significant_de_genes.all_tests.csv # Significant DE genes only across all comparisons
│   ├── fea
│   └── pseudobulk
├── annotation
├── integration
├── qc
├── report                                   # UPDATED: Interactive HTML report
└── run_info
```

Inspect `cohort.p_val_dist.png` first — a roughly uniform p-value distribution (with a spike near 0) indicates a well-powered test. Strong inflation or deflation suggests a problem with sample size or grouping.

### Functional enrichment analysis

Two complementary approaches identify enriched biological pathways from the DE results. **ORA** (Over-Representation Analysis) tests whether significantly DE genes overlap with known gene sets more than expected by chance. **GSEA** ranks all tested genes by expression and evaluates whether pathway members cluster at the top or bottom of the ranking, capturing directional pathway activity even when individual genes fall below significance thresholds. Redundant pathways are collapsed using affinity propagation (AP) clustering, and a BH-corrected FDR is applied across all comparisons.

Outputs:

```
results/
├── analysis
│   ├── differential_expression
│   ├── fea                                  # NEW
│   │   ├── gsea
│   │   │   ├── cohort.gsea.full.csv         # All comparisons combined with BH-FDR recalculated across all tests
│   │   │   ├── cohort.gsea.full.reduced.csv  # One representative gene set per AP cluster (for summary plots)
│   │   │   ├── plots
│   │   │   │   └── cohort.gsea.{comparison}.png  # NES bar chart per comparison (orange = upregulated, blue = downregulated)
│   │   │   └── male_B-cell_vs_female_B-cell      # Per-comparison WebGestaltR GSEA outputs (one folder per comparison)
│   │   │       ├── cohort.gsea.{comparison}.csv  # GSEA results: NES, p-value, core enrichment genes, AP cluster assignment
│   │   │       └── Project_gsea                  # Raw WebGestaltR outputs (enrichment table, AP clusters, HTML report)
│   │   └── ora
│   │       ├── cohort.ora.full.csv          # All comparisons combined with BH-FDR recalculated across all tests
│   │       ├── cohort.ora.full.reduced.csv   # One representative pathway per AP cluster (for summary plots)
│   │       ├── plots
│   │       │   └── cohort.ora.{comparison}.png   # log2(enrichmentRatio) bar chart per comparison
│   │       ├── male_B-cell_vs_female_B-cell      # Per-comparison WebGestaltR ORA outputs (one folder per comparison)
│   │       ├── male_Monocyte_vs_female_Monocyte
│   │       ├── male_NK-cell_vs_female_NK-cell
│   │       └── male_T-cells_vs_female_T-cells
│   │           ├── cohort.ora.{comparison}.csv   # ORA results: enrichment ratio, overlap genes, AP cluster assignment
│   │           └── Project_ora                   # Raw WebGestaltR outputs (enrichment table, AP clusters, HTML report)
│   └── pseudobulk
├── annotation
├── integration
├── qc
├── report                                   # UPDATED: Interactive HTML report
└── run_info
```

## Developer docs

For help with testing, interactive development with Singularity containers, and instructions for keeping the schema and this README in sync when parameters change:

→ [docs/developer-instructions.md](docs/developer-instructions.md)

## Component tools

The pipeline is built around the [Seurat](https://satijalab.org/seurat/) - an R framework for processing and analysing scRNAseq data. Additionally, the following R-based tools are used:

- [DoubletFinder]() for doublet detection
- [SingleR]() and [celldex]() for cell type annotation
- [DESeq2]() for the underlying differential expression analysis
- [WebGestaltR]() for both gene set enrichment analysis and overrepresentation analysis

## Additional notes

## Help / FAQ / Troubleshooting

## License(s)

## Acknowledgements/citations/credits
