process GSEA {
    publishDir "${params.outdir}/analysis/fea/gsea", mode: 'copy'
    // ext input_size: { new InputFileSizes(de_rds_path) + new InputFileSizes(gsea_db_file) }
    ext input_size: { de_rds_path.size() + ( gsea_db_file ? gsea_db_file.size() : 0 ) }
    memory { 8.GB + task.ext.input_size.B * 10 }

    input:
    tuple val(cohort_name), path(de_rds_path), val(ref_group), val(test_group), val(species)
    path gsea_db_file

    output:
    tuple val(cohort_name), path("${cohort_name}.gsea.${test_group}.${ref_group}.Rds"), emit: gsea_rds, optional: true
    tuple val(cohort_name), path("${cohort_name}.gsea.${test_group}.${ref_group}.csv"), val(ref_group), val(test_group), emit: gsea_csv, optional: true
    tuple val(cohort_name), path("${test_group}_vs_${ref_group}"), val(ref_group), val(test_group), emit: gsea_dir, optional: true


    script:
    def gsea_db_param = !!gsea_db_file ? gsea_db_file : ''
    """
    gsea.R "${cohort_name}" "${de_rds_path}" "${ref_group}" "${test_group}" "${species}" "${gsea_db_param}"
    """
}