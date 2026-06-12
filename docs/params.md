# scrnavigator-nf pipeline parameters

Pipeline parameter schema for the scrnavigator-nf Nextflow workflow.

## Input/output options

Parameters that control the pipeline inputs and outputs.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `input` | Path to the input samplesheet CSV file. | `string` |  | True |  |
| `outdir` | Directory to place output files in. | `string` | results |  |  |

## Run options

Parameters that control the pipeline execution and flow.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `qc_only` | Only run QC analysis. | `boolean` |  |  |  |
| `no_analysis` | Only run up to integration, don't perform pseudobulking, DE, or FEA. | `boolean` |  |  |  |
| `cohort_id` | Name for the integrated cohort. Default: 'cohort'. | `string` | cohort |  |  |
| `resolutions` | Comma-separated list of clustering resolutions (floats or integers) to use for QC. Defaults to values from 0.6 to 2.4 in steps of 0.2. | `string` | 0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4 |  |  |
| `integrated_resolution` | Final clustering resolution (float or integer) for the integrated dataset. | `number` | 1 |  |  |
| `cluster_method` | Method to use for clustering. Either 'louvain' or 'leiden' (default: 'leiden'). (accepted: `louvain`\|`leiden`) | `string` | leiden |  |  |
| `species` | The species that the samples belong to. | `string` |  | True |  |
| `no_mt` | Specifies that annotation of mitochondrial gene percentages should be skipped. <details><summary>Help</summary><small>Specify `no_mt` when running nuclei samples.</small></details>| `boolean` |  |  |  |
| `mt_gene_list` | Path to a text file containing known mitochondrial genes in the species. Not required when --species is either 'human' or 'mouse'. | `string` |  |  |  |
| `s_genes` | Path to a text file containing S gene IDs for cell cycle annotation. Not required when --species is 'human'. | `string` |  |  |  |
| `g2m_genes` | Path to a text file containing G2M gene IDs for cell cycle annotation. Not required when --species is 'human'. | `string` |  |  |  |
| `annotation_db` | Path to an RDS file containing a suitable reference dataset for cell type annotation with SingleR. Must contain a field called 'label' with the cell type labels. When --species is 'human', the Human Primary Cell Atlas (HPCA) data from the celldex package will be used by default. As per the SingleR documentaiton, this must be 'a numeric matrix of single-cell expression values where rows are genes and columns are cells. Alternatively, a SummarizedExperiment object containing such a matrix' | `string` |  |  |  |
| `min_cells_for_annotation` | The minimum number of cells required to keep a cell type annotation. | `integer` | 10 |  |  |
| `custom_marker_genes` | Path to a CSV file containing custom marker gene sets for annotation of cell types. Must contain two columns: 'cell_type' and 'gene_ids', where 'gene_ids' is a semicolon-delimited list of gene IDs. | `string` |  |  |  |
| `custom_annotation_mad_threshold` | Median absolute difference (MAD) threshold to use when determining cell type based on custom gene programs. | `number` | 1 |  |  |
| `cluster_annotation` | The annotation label to use for annotating clusters. If not supplied, defaults to the first of the following that exists: 'custom_cell_type', 'custom_cell_type.max_score', 'SingleR.annotation', 'SingleR.hpca_main', 'SingleR.mouse_main', 'SingleR.hpca_fine', 'SingleR.mouse_fine', 'Phase'. | `string` |  |  |  |
| `cell_type_proportion_threshold` | Proportion of cells in a cluster that must be of the same cell type before calling the cell type for the cluster. | `number` | 0.67 |  |  |
| `manual_cluster_annotations` | Path to a CSV file to manually set cluster annotations. Must have two columns: 'cluster' and 'cell_type'. | `string` |  |  |  |
| `pseudo_groups` | A comma-separated list of metadata variables to be used for grouping cells for pseudobulking and differential experession analysis. Defaults to the 'cluster_annotation' field created by the cluster annotation module. | `string` | cluster_annotation |  |  |
| `comparisons` | A path to a CSV file containing comparisons to be used for differential expression analysis. | `string` |  |  |  |
| `p_value_threshold` | P-value threshold for differential expression analysis. Default is 0.05 | `number` | 0.05 |  |  |
| `fc_threshold` | Fold-change threshold for differential expression analysis. Set to '0' to disable fold-change thresholding (default). | `number` | 0 |  |  |
| `ora_db` | Path to a .gmt file containing an enrichment database for use with over representation analysis with WebGestaltR. Must contain three tab-delimited columns: category ID, external link, and gene ID. | `string` |  |  |  |
| `gsea_db` | Path to a .gmt file containing an enrichment database for use with gene set enrichment analysis with WebGestaltR. Must contain three tab-delimited columns: category ID, external link, and gene ID. | `string` |  |  |  |
| `ens_db` | Path to an sqlite ensembldb for the species. <details><summary>Help</summary><small>Leave empty to use the built in ensembld e.g. humans: `EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86`; mice: `EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79`</small></details>| `string` |  |  |  |
| `ens_db_version` | Version of the ensembldb <details><summary>Help</summary><small>Used to download the specified version in DOWNLOAD_ENSDB. The default is v113 for the most recent available version for human.</small></details>| `string` | v113 |  |  |
| `vars_to_regress` | Comma-separated list of metadata variables to regress out when performing scaling and normalisation with the `SCTransform` function in `Seurat`. If MT genes have been annotated, the percentage of MT genes per sample (`percent.mt`) is automatically added to this list. | `string` |  |  |  |

## Miscellaneous parameters

Miscellaneous parameters.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `help` | Prints the pipeline help message. | `boolean` |  |  | True |
| `timestamp` | Timestamp of pipeline launch for reporting. | `string` |  |  | True |
