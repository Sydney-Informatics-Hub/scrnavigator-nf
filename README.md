# scrnavigator-nf

## Description

This Nextflow pipeline provides a standardised method for processing and analysing single cell RNA sequencing (scRNAseq) data. It builds upon the [scRNAvigator notebooks](https://github.com/Sydney-Informatics-Hub/scrna-analysis) - a series of interactive Quarto notebooks that step the user through quality control, dataset integration, cell annotation, pseudobulking, differential expression and functional enrichment analysis (FEA). In this Nextflow pipeline, these same steps have been organised into a semi-automated workflow that can be easily scaled up and take advantage of parallelisation and high performance computing infrastructures. The workflow is designed to be run in an iterative way, starting with quality control and filtering, then integration and cell annotation, and finally pseudobulking, differential expression and FEA. This allows the user to inspect vital quality metrics at each step and select filtering thresholds, annotations and sample groupings required for the downstream stages.

## User guide

### Input data


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
