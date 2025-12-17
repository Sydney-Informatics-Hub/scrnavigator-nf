process ANALYSE_GSEA {
    publishDir "${params.outdir}/analysis/fea/gsea", mode: 'copy'

    input:
    tuple val(cohort_name), path(full_rds_paths), path(reduced_rds_paths)

    output:
    tuple val(cohort_name), path("${cohort_name}.gsea.full.csv"), emit: gsea_csv
    tuple val(cohort_name, path("${cohort_name}.gsea.full.reduced.csv")), emit: gsea_reduced_csv

    script:
    def all_full_rds_files = full_rds_paths.join('\n')
    def all_reduced_rds_files = reduced_rds_paths.join('\n')
    """
    echo -e "${all_full_rds_files}" > all_full_rds_files.txt
    echo -e "${all_reduced_rds_files}" > all_reduced_rds_files.txt
    analyse_gsea.R "${cohort_name}" all_full_rds_files.txt all_reduced_rds_files.txt
    """
}