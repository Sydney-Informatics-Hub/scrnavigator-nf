process CLUSTER_INTEGRATED {
    publishDir "${params.outdir}/integration/clustering", mode: 'copy'
    ext input_size: { new InputFileSizes(rds_path) }
    memory { task.ext.input_size.getSizeMB() * 2 + 1.GB }

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