process INTEGRATION {
    publishDir "${params.outdir}/integration", mode: 'copy'

    input:
    path rds_paths
    val cohort_name

    output:
    tuple val(cohort_name), path("${cohort_name}.integrated.rds"), emit: integrated_rds
    tuple val(cohort_name), path("qc_results")

    script:
    def all_rds_paths = rds_paths.collect { f -> "'${f}'" }.join(' ')
    """
    # Create parameter samplesheet
    integrate.R "${cohort_name}" ${all_rds_paths}
    """
}