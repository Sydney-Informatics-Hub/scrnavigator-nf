process REPORT {
    publishDir "${params.outdir}/report", mode: 'copy'
    // ext input_size: { new InputFileSizes(rds_path) }
    // ext input_size: { rds_path.size() }
    memory 4.GB
    container "sydneyinformaticshub/scrnavigator-nf-report"

    input:
    tuple val(cohort_name),
        val(metadata),
        path(report_templates),
        path(qc_filter_results, stageAs: "qc_filter/??"),
        path(qc_cluster_results, stageAs: "qc_cluster/??"),
        path(qc_doublet_results, stageAs: "qc_doublets/??"),
        path(integration_qc_results, stageAs: "integration_qc"),
        path(integration_cluster_results, stageAs: "integration_cluster"),
        path(annotation_cc_results, stageAs: "annotation_cc"),
        path(annotation_db_results, stageAs: "annotation_db"),
        path(annotation_custom_results, stageAs: "annotation_custom"),
        path(annotation_clusters_results, stageAs: "annotation_clusters"),
        path(analysis_pseudo_comparison_groups, stageAs: "analysis_pseudo/*"),
        path(analysis_de_results, stageAs: "analysis_de/*"),
        path(analysis_gsea_results, stageAs: "analysis_gsea/*"),
        path(analysis_ora_results, stageAs: "analysis_ora/*"),
        path(available_annotation_files, stageAs: "available_annotations/??.txt"),
        path(report_style)
    path samplesheet_csv, stageAs: "samplesheets/samplesheet.csv"
    path custom_marker_genes_csv, stageAs: "samplesheets/custom_marker_genes.csv"
    path comparisons_csv, stageAs: "samplesheets/comparisons.csv"
    path quarto_yaml
    val env_vars

    output:
    tuple val(cohort_name), path("report.${cohort_name}.html"), emit: report

    script:
    def int_res = (metadata.integration_resolution == null) ? '' : metadata.integration_resolution
    def pseudo_groups = (metadata.pseudo_groups == null) ? '' : metadata.pseudo_groups.tokenize(',').join('\n')
    def cluster_annotation = (metadata.cluster_annotation == null) ? '' : metadata.cluster_annotation
    def meta_fields = (metadata.meta_fields == null) ? '' : metadata.meta_fields.join('\n')
    def param_tsv = 'parameter\tvalue\n' + env_vars.params.collect { p, v -> "${p}\t${v}" }.join('\n') + '\n'
    def software_versions_csv = 'tool,version\n' + env_vars.software_versions.collect { p, v -> "${p},${v}" }.join('\n') + '\n'
    """
    # Save important metadata to file for reporting
    echo "${cohort_name}" > cohort_name.txt
    mkdir -p available_annotations
    echo "${int_res}" > integration_resolution.txt
    echo -e "${pseudo_groups}" > pseudo_groups.txt
    echo -e "${cluster_annotation}" > cluster_annotation.txt
    echo -e "${meta_fields}" > available_annotations/metadata.txt
    if [ -d "annotation_clusters" ]; then
        echo "cluster_annotation" > available_annotations/clusters.txt
    fi
    echo -e "${param_tsv}" > parameters.tsv
    echo -e "${software_versions_csv}" > software.csv
    echo "${env_vars.pipe_version}" > version.txt
    echo "${env_vars.nxf_version}" > nextflow.txt
    echo -e "${env_vars.cmd}" > cmd.sh

    # Create a cache dir for quarto
    mkdir -p .cache
    export XDG_CACHE_HOME="\${PWD}/.cache"

    # Render the quarto website
    quarto render index.qmd --output report.${cohort_name}.html
    """
}