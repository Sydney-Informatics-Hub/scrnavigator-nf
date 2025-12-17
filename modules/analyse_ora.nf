process ANALYSE_ORA {
    publishDir "${params.outdir}/analysis/fea/ora", mode: 'copy'

    input:
    tuple val(cohort_name), path(full_rds_paths), path(reduced_rds_paths)

    output:
    tuple val(cohort_name), path("${cohort_name}.osa.full.csv"), emit: osa_csv
    tuple val(cohort_name, path("${cohort_name}.osa.full.reduced.csv")), emit: osa_reduced_csv

    script:
    def all_full_rds_files = full_rds_paths.join('\n')
    def all_reduced_rds_files = reduced_rds_paths.join('\n')
    """
    echo -e "${all_full_rds_files}" > all_full_rds_files.txt
    echo -e "${all_reduced_rds_files}" > all_reduced_rds_files.txt
    analyse_ora.R "${cohort_name}" all_full_rds_files.txt all_reduced_rds_files.txt
    """
}