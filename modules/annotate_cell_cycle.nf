process ANNOTATE_CELL_CYCLE {
    publishDir "${params.outdir}/annotation/cell_cycle", mode: 'copy'

    input:
    tuple val(cohort_name), path(rds_path)
    val annotation_params

    output:
    tuple val(cohort_name), path("${cohort_name}.annotated.cell_cycle.rds"), emit: annotated_rds
    tuple val(cohort_name), path("available_annotations.txt"), emit: available_annotations
    tuple val(cohort_name), path("qc_results"), emit: qc_results

    script:
    def annotation_csv = 'param,value\n' +
        annotation_params.collect { param, value -> [ param, ( value == null ? '' : value ) ].join(',') }.join('\n')
    """
    # Create annotation samplesheet
    echo -e "${annotation_csv}" > annotation.csv

    annotate_cell_cycle.R "${cohort_name}" "${rds_path}" annotation.csv
    """
}