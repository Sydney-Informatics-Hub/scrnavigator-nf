process CLUSTER_INTEGRATED {
    publishDir "${params.outdir}/integration/clustering", mode: 'copy'

    input:
    tuple val(cohort_name), path(rds_path), val(meta)
    

    output:
    tuple val(cohort_name), path("${cohort_name}.integrated.clustered.rds"), val(meta), emit: integrated_rds
    tuple val(cohort_name), path("qc_results"), emit: qc_results

    script:
    def params_csv = 'param,value\n' +
        meta.collect { param, value -> [ param, ( value == null ? '' : value ) ].join(',') }.join('\n')
    """
    # Create parameter samplesheet
    echo -e "${params_csv}" > params.csv

    # Create parameter samplesheet
    cluster_integrated.R "${cohort_name}" "${rds_path}" params.csv
    """
}