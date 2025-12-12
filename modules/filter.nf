process FILTER {
    publishDir "${params.outdir}/qc/${sample}/filter", mode: 'copy'

    input:
    tuple val(sample), path(rds_path), val(meta)

    output:
    tuple val(sample), path("${sample}.qc_filtered.rds"), val(meta), emit: qc_rds
    tuple val(sample), path("qc_results"), emit: qc_results

    script:
    def params_csv = 'param,value\n' +
        meta.collect { param, value -> [ param, ( value == null ? '' : value ) ].join(',') }.join('\n')
    """
    # Create parameter samplesheet
    echo -e "${params_csv}" > params.csv

    filter.R "${sample}" "${rds_path}" params.csv
    """
}