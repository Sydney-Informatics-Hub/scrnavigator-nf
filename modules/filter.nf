process FILTER {
    publishDir "${params.outdir}/qc/${sample}/filter", mode: 'copy'
    ext input_size: { new InputFileSizes(rds_path) }
    memory { task.ext.input_size.getSizeMB() * 2 + 1.GB }

    input:
    tuple val(sample), path(rds_path), val(sample_params), path(cells_to_remove)

    output:
    tuple val(sample), path("${sample}.qc_filtered.rds"), emit: qc_rds
    tuple val(sample), path("qc_results"), emit: qc_results

    script:
    def params_csv = 'param,value\n' +
        sample_params.collect { param, value -> [ param, ( value == null ? '' : value ) ].join(',') }.join('\n')
    """
    # Create parameter samplesheet
    echo -e "${params_csv}" > params.csv

    filter.R "${sample}" "${rds_path}" params.csv "${cells_to_remove}"
    """
}