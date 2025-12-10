process INIT_QC {
    input:
    tuple val(sample), path(rds_path), val(meta)
    val resolutions
    path report_template

    output:
    tuple val(sample), path("${sample}.qc.rds"), emit: qc_rds

    script:
    def params_csv = 'param,value\n' +
        meta.collect { param, value -> [ param, value ].join(',') }.join('\n') +
        "\nresolutions,${resolutions}\n"
    """
    # Create parameter samplesheet
    echo -e "${params_csv}" > params.csv

    init_qc.R "${sample}" "${rds_path}" params.csv ${report_template}
    """
}