process PSEUDOBULK {
    publishDir "${params.outdir}/analysis/pseudobulk", mode: 'copy'
    // ext input_size: { new InputFileSizes(rds_path) }
    ext input_size: { rds_path.size() }
    memory { task.ext.input_size.B * 2 + 1.GB }

    input:
    tuple val(cohort_name), path(rds_path), val(grouping_fields)

    output:
    tuple val(cohort_name), path("${cohort_name}.pseudobulk.rds"), emit: pseudobulked_rds
    tuple val(cohort_name), path("${cohort_name}.comparison_groups.txt"), emit: comparison_groups

    script:
    """
    pseudobulk.R "${cohort_name}" "${rds_path}" "${grouping_fields}"
    """
}