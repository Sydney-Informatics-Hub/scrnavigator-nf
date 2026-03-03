process CLUSTER_INTEGRATED {
    publishDir "${params.outdir}/integration/clustering", mode: 'copy'
    // ext input_size: { new InputFileSizes(rds_path) }
    ext input_size: { rds_path.size() }
    memory { 16.GB + task.ext.input_size.B * 10 }

    input:
    tuple val(cohort_name), path(rds_path), val(all_resolutions), val(cluster_method), val(integrated_resolution)

    output:
    tuple val(cohort_name), path("${cohort_name}.integrated.clustered.rds"), emit: integrated_rds
    tuple val(cohort_name), path("qc_results"), emit: qc_results

    script:
    def params_csv = 'param,value\n' +
        "resolutions,${all_resolutions}\n" +
        "cluster_method,${cluster_method}\n" +
        "integrated_resolution,${integrated_resolution}\n"
    """
    # Create parameter samplesheet
    echo -e "${params_csv}" > params.csv

    cluster_integrated.R "${cohort_name}" "${rds_path}" params.csv
    """
}