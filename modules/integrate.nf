process INTEGRATE {
    publishDir "${params.outdir}/integration", mode: 'copy'

    input:
    path rds_paths
    val cohort_name
    val resolutions
    val integrated_resolution
    val cluster_method

    output:
    tuple val(cohort_name), path("${cohort_name}.integrated.rds"), emit: integrated_rds
    tuple val(cohort_name), path("qc_results"), emit: qc_results

    script:
    """
    # Create parameter samplesheet
    integrate.R "${cohort_name}" "${rds_paths}"
    """
}