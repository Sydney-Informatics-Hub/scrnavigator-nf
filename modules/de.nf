process DIFFERENTIAL_EXPRESSION {
    publishDir "${params.outdir}/analysis/differential_expression", mode: 'copy'
    ext input_size: { new InputFileSizes(rds_path) }
    memory { task.ext.input_size.getSizeMB() * 2 + 1.GB }

    input:
    tuple val(cohort_name), path(rds_path), val(ref_group), val(test_group)

    output:
    tuple val(cohort_name), path("${cohort_name}.de.${test_group}.${ref_group}.Rds"), emit: de_rds
    tuple val(cohort_name), path("${cohort_name}.de.${test_group}.${ref_group}.csv"), val(ref_group), val(test_group), emit: de_csv

    script:
    """
    de.R "${cohort_name}" "${rds_path}" "${ref_group}" "${test_group}"
    """
}