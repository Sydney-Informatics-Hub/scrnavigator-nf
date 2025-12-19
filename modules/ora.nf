process ORA {
    publishDir "${params.outdir}/analysis/fea/ora", mode: 'copy'

    input:
    tuple val(cohort_name), path(de_rds_path), val(ref_group), val(test_group), val(species)
    tuple path(ora_db_file)

    output:
    tuple val(cohort_name), path("${cohort_name}.ora.${test_group}.${ref_group}.Rds"), emit: ora_rds, optional: true
    tuple val(cohort_name), path("${cohort_name}.ora.${test_group}.${ref_group}.csv"), val(ref_group), val(test_group), emit: ora_csv, optional: true
    tuple val(cohort_name), path("${test_group}_vs_${ref_group}"), val(ref_group), val(test_group), emit: ora_dir, optional: true

    script:
    def ora_db_param = !!ora_db_file ? ora_db_file : ''
    """
    ora.R "${cohort_name}" "${de_rds_path}" "${ref_group}" "${test_group}" "${species}" "${ora_db_param}"
    """
}