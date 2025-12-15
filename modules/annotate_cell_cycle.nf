process ANNOTATE_CELL_CYCLE {
    publishDir "${params.outdir}/annotation/cell_cycle", mode: 'copy'

    input:
    tuple val(cohort_name), path(rds_path), val(meta)
    val species_params

    output:
    tuple val(cohort_name), path("${cohort_name}.annotated.cell_cycle.rds"), val(meta), emit: annotated_rds
    tuple val(cohort_name), path("qc_results")

    script:
    def species_csv = 'param,value\n' +
        species_params.collect { param, value -> [ param, ( value == null ? '' : value ) ].join(',') }.join('\n')
    """
    # Create species samplesheet
    echo -e "${species_csv}" > species.csv

    annotate_cell_cycle.R "${cohort_name}" "${rds_path}" species.csv
    """
}