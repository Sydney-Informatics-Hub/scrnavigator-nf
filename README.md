# scrnavigator-nf

## Description

This Nextflow pipeline provides a standardised method for processing and analysing single cell RNA sequencing (scRNAseq) data. It builds upon the [scRNAvigator notebooks](https://github.com/Sydney-Informatics-Hub/scrna-analysis) - a series of interactive Quarto notebooks that step the user through quality control, dataset integration, cell annotation, pseudobulking, differential expression and functional enrichment analysis (FEA). In this Nextflow pipeline, these same steps have been organised into a semi-automated workflow that can be easily scaled up and take advantage of parallelisation and high performance computing infrastructures. The workflow is designed to be run in an iterative way, starting with quality control and filtering, then integration and cell annotation, and finally pseudobulking, differential expression and FEA. This allows the user to inspect vital quality metrics at each step and select filtering thresholds, annotations and sample groupings required for the downstream stages.

## User guide

### Pipeline setup

```bash
git clone git@github.com:Sydney-Informatics-Hub/scrnavigator-nf.git
cd scrnavigator-nf
```

### Input data

The initial input data that this pipeline requires is one or more RDS files, each containing a Seurat object with the count matrix data for a single sample. [Seurat is an R package](https://satijalab.org/seurat/) that provides a framework for handling and processing single cell RNAseq data.

The initial data is expected to have already been pre-processed to remove background noise, i.e. ambient RNA signals. This is done automatically by 10X's `cellranger` software for their platform, and similar filtering can be achieved with tools like [`cellbender`](https://github.com/broadinstitute/CellBender) for other datasets.

The Seurat data object should at the very minimum contain the count matrix data for your sample in its `RNA` assay, along with either Ensembl gene IDs or HGNC gene symbols as the `RNA` assay's row names. The Seurat object's `meta.data` field should also be present with a row for every cell in the sample.

This pipeline was designed to work with the filtered output from the [nf-core `scrnaseq` Nextflow pipeline](https://nf-co.re/scrnaseq/). The `scrnaseq` pipeline can handle various scRNAseq datasets, including those from the 10X platform, and will generate RDS files containing the Seurat data required for this pipeline. We highly recommend using the `scrnaseq` pipeline for the initial alignment, counting and pre-processing of your data prior to using this workflow.

### Example samplesheet

## Parameters

Pipeline parameter schema for the scrnavigator-nf Nextflow workflow. Generated using `nf-core pipelines schema docs`.

### Input/output options 
Parameters that control the pipeline inputs and outputs.                                                                                        
| Parameter | Description | Type | Default | Required | Hidden |                                                                                
|-----------|-----------|-----------|-----------|-----------|-----------|                                                                       
| `input` | Path to the input samplesheet CSV file. | `string` |  | True |  |                                                                   
| `outdir` | Directory to place output files in. | `string` | results |  |  |                                                                   
### Run options                                                                                                                                  
Parameters that control the pipeline execution and flow.                                                                                        
| Parameter | Description | Type | Default | Required | Hidden |                                                                                
|-----------|-----------|-----------|-----------|-----------|-----------|                                                                       
| `qc_only` | Only run QC analysis. | `boolean` | False |  |  |                                                                                 
| `no_analysis` | Only run up to integration, don't perform pseudobulking, DE, or FEA. | `boolean` | False |  |  |                              
| `cohort_id` | Name for the integrated cohort. Default: 'cohort'. | `string` | cohort |  |  |                                                  
| `resolutions` | Comma-separated list of clustering resolutions (floats or integers) to use for QC. Defaults to values from 0.6 to 2.4 in steps of 0.2. | `string` |  |  |  |             
| `integrated_resolution` | Final clustering resolution (float or integer) for the integrated dataset. | `number` |  |  |  |                    
| `cluster_method` | Method to use for clustering. Either 'louvain' or 'leiden' (default: 'leiden'). (accepted: `louvain`\|`leiden`) | `string` 
| leiden |  |  |
| `species` | The species that the samples belong to. | `string` |  |  |  |
| `no_mt` | Specifies that annotation of mitochondrial gene percentages should be skipped. | `boolean` | False |  |  |
| `mt_gene_list` | Path to a text file containing known mitochondrial genes in the species. Not required when --species is either 'human' or 'mouse'. | `string` | None |  |  |
| `ens_db_rds` | Path to an RDS file containing an AnnotationDb object for the Ensembl database (e.g. EnsDb.Hsapiens.v86). Required for species other than 'human' or 'mouse'. | `string` | None |  |  |
| `s_genes` | Path to a text file containing S gene IDs for cell cycle annotation. Not required when --species is 'human'. | `string` | None |
| `g2m_genes` | Path to a text file containing G2M gene IDs for cell cycle annotation. Not required when --species is 'human'. | `string` | None |  |  |  
| `annotation_db` | Path to an RDS file containing a suitable reference dataset for cell type annotation with SingleR. Must contain a field called 'label' with the cell type labels. When --species is 'human', the Human Primary Cell Atlas (HPCA) data from the celldex package will be used by default. As per the SingleR documentaiton, this must be 'a numeric matrix of single-cell expression values where rows are genes and columns are cells. Alternatively, a SummarizedExperiment object containing such a matrix' | `string` | None |  |  |
| `min_cells_for_annotation` | The minimum number of cells required to keep a cell type annotation. | `integer` | 10 |  |  |
| `custom_marker_genes` | Path to a CSV file containing custom marker gene sets for annotation of cell types. Must contain two columns: 'cell_type' and 'gene_ids', where 'gene_ids' is a semicolon-delimited list of gene IDs. | `string` | None |  |  |
| `custom_annotation_mad_threshold` | Median absolute difference (MAD) threshold to use when determining cell type based on custom gene programs. | `number` | 1 |  |  |
| `cluster_annotation` | The annotation label to use for annotating clusters. If not supplied, defaults to the first of the following that exists: 'custom_cell_type', 'custom_cell_type.max_score', 'SingleR.annotation', 'SingleR.hpca_main', 'SingleR.hpca_fine', 'Phase'. | `string` | None |  |  |
| `cell_type_proportion_threshold` | Proportion of cells in a cluster that must be of the same cell type before calling the cell type for the cluster. | `number` | 0.67 |  |  |
| `manual_cluster_annotations` | Path to a CSV file to manually set cluster annotations. Must have two columns: 'cluster' and 'cell_type'. | `string` | None |  |  |
| `pseudo_groups` | A comma-separated list of metadata variables to be used for grouping cells for pseudobulking and differential experession analysis. Defaults to the 'cluster_annotation' field created by the cluster annotation module. | `string` | cluster_annotation |  |  |
| `comparisons` | A path to a CSV file containing comparisons to be used for differential expression analysis. | `string` | None |  |  |
| `p_value_threshold` | P-value threshold for differential expression analysis. Default is 0.05 | `number` | 0.05 |  |  |
| `fc_threshold` | Fold-change threshold for differential expression analysis. No threshold is applied by default. | `number` | None |  |  |
| `ora_db` | Path to a .gmt file containing an enrichment database for use with over representation analysis with WebGestaltR. Must contain three tab-delimited columns: category ID, external link, and gene ID. | `string` | None |  |  |
| `gsea_db` | Path to a .gmt file containing an enrichment database for use with gene set enrichment analysis with WebGestaltR. Must contain three tab-delimited columns: category ID, external link, and gene ID. | `string` | None |  |  |

### Miscellaneous parameters                                                                                                                     
| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|                                                                       
| `help` | Prints the pipeline help message. | `boolean` |  |  | True |

## Usage

NOTE: The pipeline is currently under active development and usage may change rapidly.

While this pipeline can be run end-to-end, the steps are recommended to be run sequentially, inspecting the data and outputs of a subworkflow to inform the parameters to use in the subsequent steps. 

We recommend reviewing the `index.html` report after each step has run, which compiles the results of each step in an interactive report for convenient EDA.

### Quality control and filtering

Run pre-processing and quality control with no cell filtering. Good first pass to understand your data for further filtering. The `--qc_only` option conducts only quality control steps to inspect and remove low-quality cells for each sample defined in the samplesheet.

```bash
nextflow run main.nf \
    --input path/to/samplesheet.csv \
    --outdir results \
    --species human \
    --qc_only \
    --cluster_method louvain \
    -profile gadi \
    --gadi_account er01 \
    --gadi_storage scratch/er01+gdata/er01+gdata/if89 \
    -c /scratch/er01/PIPE-6945-scrna/test.scrnavigator.opt.config
```

Example minimal samplesheet:

```console
sample,rds,sex
Donor_1,/g/data/er01/test-data/single-cell/human/rds/Donor1_filtered_matrix.seurat.rds
Donor_2,/g/data/er01/test-data/single-cell/human/rds/Donor2_filtered_matrix.seurat.rds
Donor_3,/g/data/er01/test-data/single-cell/human/rds/Donor3_filtered_matrix.seurat.rds
```

Output:

```
results                
├── qc # These contain the intermediate QC results for each sample that are used to generate the report (below)
│   ├── Donor_1
│   │   ├── cluster # Clustering plots, barcode cell assignment across resolutions  
│   │   ├── filter  # Low-quality cell filtering metrics
│   │   └── preprocess
│   ├── Donor_2
│   │   └── ...
├── report
│   └── report # Interactive HTML reports to explore QC and clustering results
└── run_info
```

From inspecting the report from the previous output, identify filtering thresholds for low-quality cells based on the:

* Number of RNA counts (3k - 12k)
* Features (> 2000)
* Percentage of mitochondrial genes mapped (< 10%)

Update the samplesheet to include these thresholds. Thresholds with no upper or lower bound are left intentionally blank:

```console
sample,rds,sex,min_ncount,max_ncount,min_nfeature,max_nfeature,min_mt_pct,max_mt_pct
Donor_1,/g/data/er01/test-data/single-cell/human/rds/Donor1_filtered_matrix.seurat.rds,M,3000,12000,2000,20000,,10,1.2
Donor_2,/g/data/er01/test-data/single-cell/human/rds/Donor2_filtered_matrix.seurat.rds,M,3000,12000,2000,20000,,10,1
Donor_3,/g/data/er01/test-data/single-cell/human/rds/Donor3_filtered_matrix.seurat.rds,F,3000,12000,2000,20000,,10,1.2
```

Re-run QC with the amended samplesheet using the same command and `-resume`:

```bash
nextflow run main.nf \
    --input path/to/samplesheet.csv \
    --outdir results \
    --species human \
    --qc_only \
    --cluster_method louvain \
    -profile gadi \
    --gadi_account er01 \
    --gadi_storage scratch/er01+gdata/er01+gdata/if89 \
    -c /scratch/er01/PIPE-6945-scrna/test.scrnavigator.opt.config \
    -resume
```

### Doublet detection

To conduct the next step after doublet detection, resolutions need to be selected and defined in the samplesheet for each sample. For example, add a new column called `res`:

```console
sample,rds,sex,min_ncount,max_ncount,min_nfeature,max_nfeature,min_mt_pct,max_mt_pct,res
Donor_1,/g/data/er01/test-data/single-cell/human/rds/Donor1_filtered_matrix.seurat.rds,M,3000,12000,2000,20000,,10,1.2
Donor_2,/g/data/er01/test-data/single-cell/human/rds/Donor2_filtered_matrix.seurat.rds,M,3000,12000,2000,20000,,10,1
Donor_3,/g/data/er01/test-data/single-cell/human/rds/Donor3_filtered_matrix.seurat.rds,F,3000,12000,2000,20000,,10,1.2
Donor_4,/g/data/er01/test-data/single-cell/human/rds/Donor4_filtered_matrix.seurat.rds,F,3000,12000,2000,20000,,10,0.6
```

To proceed with doublet detection after conducting QC, replace the `--qc_only` flag with `--no_analysis`:

```bash
nextflow run main.nf \
    --input path/to/samplesheet.csv \
    --outdir results \
    --species human \
    --no_analysis \
    --cluster_method louvain \
    -profile gadi \
    --gadi_account er01 \
    --gadi_storage scratch/er01+gdata/er01+gdata/if89 \
    -c /scratch/er01/PIPE-6945-scrna/test.scrnavigator.opt.config \
    -resume
```

This will use the cached results from the `QUALITY_CONTROL` steps and resume with the doublet detection step:

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

Output:

```
results                
├── integration # NEW: results of integrating all samples into a cohort RDS
│   ├── clustering
│   ├── cohort.integrated.rds
│   ├── gene_symbols.Rds
│   └── qc_results
├── qc  
│   ├── Donor_1        
│   │   ├── cluster
│   │   ├── doublets # NEW: Results and plots of doublet detection
│   │   ├── filter
│   │   └── preprocess 
│   ├── Donor_2
│   │   └── ...
├── report     
│   └── report # UPDATED: Interactive HTML reports to explore QC and clustering results
└── run_info
```

```
doublets
├── Donor_1.doublets_detected.rds # File which identifies putative doublets, without any filtering
├── Donor_1.doublets_removed.sct_clustered.rds # File with doublets removed, re-normalised and clustered.
└── qc_results # Table and graph summaries of doublet detection, used in the interactive report.
    ├── cluster_cells
    ├── Donor_1.clustree.png
    ├── Donor_1.doublet_proportions_per_cluster.1.png
    ├── Donor_1.doublets_per_cluster.1.png
    ├── Donor_1.doublet_summary.csv
    ├── Donor_1.doublet_umap.1.png
    ├── Donor_1.feature_count_plot.cluster.0.6.png   
    ├── Donor_1.feature_count_plot.cluster.0.8.png   
    ├── Donor_1.feature_count_plot.cluster.1.2.png   
    ├── Donor_1.feature_count_plot.cluster.1.4.png   
    ├── Donor_1.feature_count_plot.cluster.1.6.png   
    ├── Donor_1.feature_count_plot.cluster.1.8.png   
    ├── Donor_1.feature_count_plot.cluster.1.png     
    ├── Donor_1.feature_count_plot.cluster.2.2.png   
    ├── Donor_1.feature_count_plot.cluster.2.4.png   
    └── Donor_1.feature_count_plot.cluster.2.png     
```

### Dataset integration

Integration is conducted as part of the previous step with the `--no_analysis` option.

Key outputs:

```
results/
├── integration
│   ├── clustering
│   │   ├── cohort.integrated.clustered.rds # Integrated Seurat object with cluster assignments at multiple resolutions
│   │   └── qc_results                      # Clustree, UMAP, and per-cluster feature/count plots
│   ├── cohort.integrated.rds               # Integrated Seurat object (integrated, UMAP computed, pre-clustering)
│   ├── gene_symbols.Rds                    # Dataframe mapping gene symbols to Ensembl IDs for the cohort
│   └── qc_results
│       ├── cohort.umap.integrated.png      # UMAP coloured by sample after integration (batch correction applied)
│       └── cohort.umap.merged.png          # UMAP coloured by sample before integration (merged only)
└── report
    └── report
        ├── index.html                      # UPDATED: Interactive report
        └── ...
```

### Cell annotation

```bash
nextflow run main.nf \
    --input path/to/samplesheet.csv \
    --outdir results \
    --species human \
    --cluster_method louvain \
    -profile gadi \
    --gadi_account er01 \
    --gadi_storage scratch/er01+gdata/er01+gdata/if89 \
    -c /scratch/er01/PIPE-6945-scrna/test.scrnavigator.opt.config \
    -resume
```
### Pseudobulking



### Differential expression analysis



### Funciontal enrichment analysis



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
