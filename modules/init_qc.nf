process INIT_QC {
    publishDir "${params.outdir}/qc/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(rds_path), val(meta)
    val resolutions
    val cluster_method

    output:
    tuple val(sample), path("${sample}.qc.rds"), emit: qc_rds
    tuple val(sample), path("qc_results"), emit: qc_results

    script:
    def params_csv = 'param,value\n' +
        meta.collect { param, value -> [ param, ( value == null ? '' : value ) ].join(',') }.join('\n') +
        "\nresolutions,${resolutions}" +
        "\ncluster_method,${cluster_method}\n"
    """
    # Create parameter samplesheet
    echo -e "${params_csv}" > params.csv

    init_qc.R "${sample}" "${rds_path}" params.csv
    """
}