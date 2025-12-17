process PSEUDOBULK {
    publishDir "${params.outdir}/analysis/pseudobulk", mode: 'copy'

    input:
    tuple val(cohort_name), path(rds_path), val(grouping_fields)

    output:
    tuple val(cohort_name), path("${cohort_name}.pseudobulk.rds"), emit: pseudobulked_rds

    script:
    """
    pseudobulk.R "${cohort_name}" "${rds_path}" "${grouping_fields}"
    """
}