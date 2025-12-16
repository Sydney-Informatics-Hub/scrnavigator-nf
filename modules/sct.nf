process SCTRANSFORM {
    publishDir "${params.outdir}/qc/${sample}/cluster", mode: 'copy'

    input:
    tuple val(sample), path(rds_path), val(default_res), val(all_resolutions), val(cluster_method)

    output:
    tuple val(sample), path("${sample}.sct_clustered.rds"), emit: sct_rds
    tuple val(sample), path("qc_results"), emit: qc_results

    script:
    def res = default_res == null ? '' : default_res
    def params_csv = 'param,value\n' +
        "res,${res}\n" +
        "resolutions,${all_resolutions}\n" +
        "cluster_method,${cluster_method}\n"
    """
    # Create parameter samplesheet
    echo -e "${params_csv}" > params.csv

    sct.R "${sample}" "${rds_path}" params.csv
    """
}