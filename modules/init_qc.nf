process INIT_QC {
    input:
    tuple val(sample), path(rds_path), val(qc_param_map)
    val resolutions
    path report_template

    output:
    tuple val(sample), path("${sample}.qc.rds"), emit: qc_rds

    script:
    def min_ncount = qc_param_map.min_ncount
    def max_ncount = qc_param_map.max_ncount
    def min_nfeature = qc_param_map.min_nfeature
    def max_nfeature = qc_param_map.max_nfeature
    def min_mt_pct = qc_param_map.min_mt_pct
    def max_mt_pct = qc_param_map.max_mt_pct
    def final_resolution = qc_param_map.res
    def clusters_to_remove = qc_param_map.clusters_to_remove
    """
    # Create parameter samplesheet
    echo "param,value" > params.csv
    echo "min_ncount,${min_ncount}" >> params.csv
    echo "max_ncount,${max_ncount}" >> params.csv
    echo "min_nfeature,${min_nfeature}" >> params.csv
    echo "max_nfeature,${max_nfeature}" >> params.csv
    echo "min_mt_pct,${min_mt_pct}" >> params.csv
    echo "max_mt_pct,${max_mt_pct}" >> params.csv
    echo "resolutions,${resolutions}" >> params.csv
    echo "final_resolution,${final_resolution}" >> params.csv
    echo "clusters_to_remove,${clusters_to_remove}" >> params.csv

    init_qc.R "${sample}" "${rds_path}" params.csv ${report_template}
    """
}