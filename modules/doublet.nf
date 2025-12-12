process DETECT_DOUBLETS {
    publishDir "${params.outdir}/qc/${sample}/doublets", mode: 'copy'

    input:
    tuple val(sample), path(rds_path, stageAs: "input/*"), val(meta)
    val resolutions
    val cluster_method

    output:
    tuple val(sample), path("${sample}.doublets_removed.sct_clustered.rds"), val(meta), emit: doublets_removed_rds
    tuple val(sample), path("${sample}.doublets_detected.rds"), val(meta), emit: doublets_marked_rds
    tuple val(sample), path("qc_results"), emit: qc_results

    script:
    def params_csv = 'param,value\n' +
        meta.collect { param, value -> [ param, ( value == null ? '' : value ) ].join(',') }.join('\n') +
        "\nresolutions,${resolutions}" +
        "\ncluster_method,${cluster_method}\n"
    """
    # Create parameter samplesheet
    echo -e "${params_csv}" > params.csv

    doublet.R "${sample}" "${rds_path}" params.csv

    # Once doublets are removed, re-run SCTransform and clustering
    sct.R "${sample}" "${sample}.doublets_removed.rds" params.csv
    mv "${sample}.sct_clustered.rds" "${sample}.doublets_removed.sct_clustered.rds"
    """
}