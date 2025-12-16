process ANNOTATE_CLUSTERS {
    publishDir "${params.outdir}/annotation/clusters", mode: 'copy'

    input:
    tuple val(cohort_name), path(rds_path), path(manual_cluster_annotations)
    val annotation_params

    output:
    tuple val(cohort_name), path("${cohort_name}.annotated.clusters.rds"), emit: annotated_rds
    tuple val(cohort_name), path("qc_results"), emit: qc_results
    tuple val(cohort_name), path("${cohort_name}.cluster_cell_type_consensus.csv"), emit: cell_type_assignments

    script:
    def annotation_csv = 'param,value\n' +
        annotation_params.collect { param, value -> [ param, ( value == null ? '' : value ) ].join(',') }.join('\n')
    def manual_annotations_file = manual_cluster_annotations != null ? manual_cluster_annotations : ''
    """
    # Create annotation samplesheet
    echo -e "${annotation_csv}" > annotation.csv

    annotate_clusters.R "${cohort_name}" "${rds_path}" annotation.csv ${manual_annotations_file}
    """
}