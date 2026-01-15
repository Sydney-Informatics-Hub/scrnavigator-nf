process ANALYSE_DE {
    publishDir "${params.outdir}/analysis/differential_expression", mode: 'copy'
    // ext input_size: { new InputFileSizes(de_rds_paths) }
    ext input_size: { de_rds_paths instanceof Collection ? de_rds_paths.inject(0) { sum, f -> sum + f.size() } : de_rds_paths.size() }
    memory { task.ext.input_size.B * 2 + 1.GB }

    input:
    tuple val(cohort_name), path(de_rds_paths), path(gene_symbols), val(p_thresh), val(fc_thresh)

    output:
    tuple val(cohort_name), path("${cohort_name}.de.full.Rds"), emit: de_rds
    tuple val(cohort_name), path("${cohort_name}.de.full.csv"), emit: de_csv
    tuple val(cohort_name), path("results")

    script:
    def fc_thresh_param = fc_thresh == null ? '' : fc_thresh
    def de_params_csv = 'param,value\n' +
        "p_val_cutoff,${p_thresh}\n" +
        "fc_cutoff,${fc_thresh_param}\n"
    def all_rds_paths = de_rds_paths.collect { f -> "'${f}'" }.join(' ')
    """
    # Create parameter samplesheet
    echo -e "${de_params_csv}" > params.csv

    analyse_de.R "${cohort_name}" "${gene_symbols}" params.csv ${all_rds_paths}
    """
}