process ANNOTATE_DATABASE {
    publishDir "${params.outdir}/annotation/db", mode: 'copy'

    input:
    tuple val(cohort_name), path(rds_path), val(meta)
    val species_params

    output:
    tuple val(cohort_name), path("${cohort_name}.annotated.database.rds"), val(meta), emit: annotated_rds
    tuple val(cohort_name), path("qc_results")

    script:
    def species_csv = 'param,value\n' +
        species_params.collect { param, value -> [ param, ( value == null ? '' : value ) ].join(',') }.join('\n')
    """
    # Create species samplesheet
    echo -e "${species_csv}" > species.csv

    annotate_db.R "${cohort_name}" "${rds_path}" species.csv
    """
}