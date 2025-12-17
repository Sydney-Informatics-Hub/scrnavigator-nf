process ORA {
    publishDir "${params.outdir}/analysis/fea/ora", mode: 'copy'

    input:
    tuple val(cohort_name), path(de_rds_path), val(ref_group), val(test_group)

    output:
    tuple val(cohort_name), path("${cohort_name}.ora.${test_group}.${ref_group}.Rds"), emit: ora_rds, optional: true
    tuple val(cohort_name), path("${cohort_name}.ora.${test_group}.${ref_group}.csv"), val(ref_group), val(test_group), emit: ora_csv, optional: true
    tuple val(cohort_name), path("${cohort_name}.ora.${test_group}.${ref_group}.reduced.Rds"), emit: ora_reduced_rds, optional: true
    tuple val(cohort_name), path("${cohort_name}.ora.${test_group}.${ref_group}.reduced.csv"), val(ref_group), val(test_group), emit: ora_reduced_csv, optional: true

    script:
    """
    ora.R "${cohort_name}" "${de_rds_path}" "${ref_group}" "${test_group}"
    """
}