process ANNOTATE_CUSTOM {
    publishDir "${params.outdir}/annotation/custom", mode: 'copy'
    // ext input_size: { new InputFileSizes(rds_path) + new InputFileSizes(ens_db_rds) }
    ext input_size: { rds_path.size() + ( ens_db_rds ? ens_db_rds.size() : 0 ) }
    memory { 8.GB + task.ext.input_size.B * 10 }

    input:
    tuple val(cohort_name), path(rds_path), val(species), path(ens_db_rds), path(custom_marker_genes), val(custom_annotation_mad_threshold)

    output:
    tuple val(cohort_name), path("${cohort_name}.annotated.custom.rds"), emit: annotated_rds
    tuple val(cohort_name), path("available_annotations.txt"), emit: available_annotations
    tuple val(cohort_name), path("qc_results"), emit: qc_results

    script:
    def annotation_csv = 'param,value\n' +
        "species,${species}\n" +
        "ens_db_rds,${ens_db_rds}\n" +
        "custom_marker_genes,${custom_marker_genes}\n" +
        "custom_annotation_mad_threshold,${custom_annotation_mad_threshold}\n"
    """
    # Create annotation samplesheet
    echo -e "${annotation_csv}" > annotation.csv

    annotate_custom.R "${cohort_name}" "${rds_path}" annotation.csv
    """
}