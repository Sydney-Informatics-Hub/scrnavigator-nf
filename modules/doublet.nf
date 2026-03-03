process DETECT_DOUBLETS {
    publishDir "${params.outdir}/qc/${sample}/doublets", mode: 'copy'
    // ext input_size: { new InputFileSizes(rds_path) }
    ext input_size: { rds_path.size() }
    memory { 16.GB + task.ext.input_size.B * 10 }

    input:
    tuple val(sample), path(rds_path, stageAs: "input/*"), val(multiplet_rate), val(default_res), val(all_resolutions), val(cluster_method)

    output:
    tuple val(sample), path("${sample}.doublets_removed.sct_clustered.rds"), emit: doublets_removed_rds
    tuple val(sample), path("${sample}.doublets_detected.rds"), emit: doublets_marked_rds
    tuple val(sample), path("qc_results"), emit: qc_results

    script:
    def mr = multiplet_rate == null ? '' : multiplet_rate
    def res = default_res == null ? '' : default_res
    def params_csv = 'param,value\n' +
        "multiplet_rate,${mr}\n" +
        "res,${res}\n" +
        "resolutions,${all_resolutions}\n" +
        "cluster_method,${cluster_method}\n"
    """
    # Create parameter samplesheet
    echo -e "${params_csv}" > params.csv

    doublet.R "${sample}" "${rds_path}" params.csv

    # Once doublets are removed, re-run SCTransform and clustering
    sct.R "${sample}" "${sample}.doublets_removed.rds" params.csv
    mv "${sample}.sct_clustered.rds" "${sample}.doublets_removed.sct_clustered.rds"
    """
}