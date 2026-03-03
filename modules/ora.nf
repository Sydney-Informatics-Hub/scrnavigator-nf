process ORA {
    publishDir "${params.outdir}/analysis/fea/ora", mode: 'copy'
    // ext input_size: { new InputFileSizes(de_rds_path) + new InputFileSizes(ora_db_file) }
    ext input_size: { de_rds_path.size() + ( ora_db_file ? ora_db_file.size() : 0 ) }
    memory { 8.GB + task.ext.input_size.B * 10 }

    input:
    tuple val(cohort_name), path(de_rds_path), val(ref_group), val(test_group), val(species)
    path ora_db_file

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