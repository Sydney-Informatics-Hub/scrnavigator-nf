process ANNOTATE_CLUSTERS {
    publishDir "${params.outdir}/annotation/clusters", mode: 'copy'

    input:
    tuple val(cohort_name), path(rds_path), val(species), val(cluster_annotation), val(cell_type_proportion_threshold), path(manual_cluster_annotations)

    output:
    tuple val(cohort_name), path("${cohort_name}.annotated.clusters.rds"), emit: annotated_rds
    tuple val(cohort_name), path("qc_results"), emit: qc_results
    tuple val(cohort_name), path("${cohort_name}.cluster_cell_type_consensus.csv"), emit: cell_type_assignments

    script:
    def annotation_csv = 'param,value\n' +
        "species,${species}\n" +
        "cluster_annotation,${cluster_annotation}\n" +
        "cell_type_proportion_threshold,${cell_type_proportion_threshold}\n" +
        "manual_cluster_annotations,${manual_cluster_annotations}\n"
    """
    # Create annotation samplesheet
    echo -e "${annotation_csv}" > annotation.csv

    annotate_clusters.R "${cohort_name}" "${rds_path}" annotation.csv
    """
}