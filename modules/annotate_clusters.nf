process ANNOTATE_CLUSTERS {
    publishDir "${params.outdir}/annotation/clusters", mode: 'copy'

    input:
    tuple val(cohort_name), path(rds_path), val(species), val(cluster_annotation), val(cell_type_proportion_threshold), path(manual_cluster_annotations)

    output:
    tuple val(cohort_name), path("${cohort_name}.annotated.clusters.rds"), emit: annotated_rds
    tuple val(cohort_name), path("qc_results"), emit: qc_results
    tuple val(cohort_name), path("${cohort_name}.cluster_cell_type_assignments.csv"), emit: cell_type_assignments

    script:
    def cluster_ann = cluster_annotation == null ? '' : cluster_annotation
    def man_clusters = manual_cluster_annotations == null ? '' : manual_cluster_annotations
    def annotation_csv = 'param,value\n' +
        "species,${species}\n" +
        "cluster_annotation,${cluster_ann}\n" +
        "cell_type_proportion_threshold,${cell_type_proportion_threshold}\n" +
        "manual_cluster_annotations,${man_clusters}\n"
    """
    # Create annotation samplesheet
    echo -e "${annotation_csv}" > annotation.csv

    annotate_clusters.R "${cohort_name}" "${rds_path}" annotation.csv
    """
}