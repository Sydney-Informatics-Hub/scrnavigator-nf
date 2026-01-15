process ANALYSE_ORA {
    publishDir "${params.outdir}/analysis/fea/ora", mode: 'copy'
    // ext input_size: { new InputFileSizes(rds_paths) }
    ext input_size: { rds_paths instanceof Collection ? rds_paths.inject(0) { sum, f -> sum + f.size() } : rds_paths.size() }
    memory { task.ext.input_size.B * 2 + 1.GB }

    input:
    tuple val(cohort_name), path(rds_paths)

    output:
    tuple val(cohort_name), path("${cohort_name}.ora.full.csv"), emit: ora_csv
    tuple val(cohort_name), path("${cohort_name}.ora.full.reduced.csv"), emit: ora_reduced_csv
    tuple val(cohort_name), path("plots"), emit: ora_plots

    script:
    def all_rds_paths = rds_paths.collect { f -> "'${f}'" }.join(' ')
    """
    analyse_ora.R "${cohort_name}" ${all_rds_paths}
    """
}