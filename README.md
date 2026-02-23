# scrnavigator-nf

## Description

This Nextflow pipeline provides a standardised method for processing and analysing single cell RNA sequencing (scRNAseq) data. It builds upon the [scRNAvigator notebooks](https://github.com/Sydney-Informatics-Hub/scrna-analysis) - a series of interactive Quarto notebooks that step the user through quality control, dataset integration, cell annotation, pseudobulking, differential expression and functional enrichment analysis (FEA). In this Nextflow pipeline, these same steps have been organised into a semi-automated workflow that can be easily scaled up and take advantage of parallelisation and high performance computing infrastructures. The workflow is designed to be run in an iterative way, starting with quality control and filtering, then integration and cell annotation, and finally pseudobulking, differential expression and FEA. This allows the user to inspect vital quality metrics at each step and select filtering thresholds, annotations and sample groupings required for the downstream stages.

## User guide

### Input data

The initial input data that this pipeline requires is one or more RDS files, each containing a Seurat object with the count matrix data for a single sample. [Seurat is an R package](https://satijalab.org/seurat/) that provides a framework for handling and processing single cell RNAseq data.

The initial data is expected to have already been pre-processed to remove background noise, i.e. ambient RNA signals. This is done automatically by 10X's `cellranger` software for their platform, and similar filtering can be achieved with tools like [`cellbender`](https://github.com/broadinstitute/CellBender) for other datasets.

The Seurat data object should at the very minimum contain the count matrix data for your sample in its `RNA` assay, along with either Ensembl gene IDs or HGNC gene symbols as the `RNA` assay's row names. The Seurat object's `meta.data` field should also be present with a row for every cell in the sample.

This pipeline was designed to work with the filtered output from the [nf-core `scrnaseq` Nextflow pipeline](https://nf-co.re/scrnaseq/). The `scrnaseq` pipeline can handle various scRNAseq datasets, including those from the 10X platform, and will generate RDS files containing the Seurat data required for this pipeline. We highly recommend using the `scrnaseq` pipeline for the initial alignment, counting and pre-processing of your data prior to using this workflow.

### Pipeline setup



### Quality control and filtering



### Doublet detection



### Dataset integration



### Cell annotation



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
