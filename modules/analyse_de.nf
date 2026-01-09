process ANALYSE_DE {
    publishDir "${params.outdir}/analysis/differential_expression", mode: 'copy'

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
    """
    # Create parameter samplesheet
    echo -e "${de_params_csv}" > params.csv

    analyse_de.R "${cohort_name}" "${gene_symbols}" params.csv "${de_rds_paths}"
    """
}