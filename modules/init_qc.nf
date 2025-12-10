process INIT_QC {
    input:
    tuple val(sample), path(rds_path), val(meta)
    val resolutions
    val cluster_method

    output:
    tuple val(sample), path("${sample}.qc.rds"), emit: qc_rds

    script:
    def params_csv = 'param,value\n' +
        meta.collect { param, value -> [ param, value ].join(',') }.join('\n') +
        "\nresolutions,${resolutions}" +
        "\ncluster_method,${cluster_method}\n"
    """
    # Create parameter samplesheet
    echo -e "${params_csv}" > params.csv

    init_qc.R "${sample}" "${rds_path}" params.csv
    """
}