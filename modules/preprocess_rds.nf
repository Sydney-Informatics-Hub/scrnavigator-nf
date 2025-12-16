process PREPROCESS_RDS {
    publishDir "${params.outdir}/qc/${sample}/preprocess", mode: 'copy'

    input:
    tuple val(sample), path(rds_path), val(meta)
    val annotation_params

    output:
    tuple val(sample), path("${sample}.preprocessed.rds"), emit: preprocessed_rds

    script:
    def meta_csv = 'field,value\n' +
        meta.collect { field, value -> [ field, ( value == null ? '' : value ) ].join(',') }.join('\n')
    def annotation_csv = 'param,value\n' +
        annotation_params.collect { param, value -> [ param, ( value == null ? '' : value ) ].join(',') }.join('\n')
    """
    # Create parameter samplesheet
    echo -e "${meta_csv}" > meta.csv

    # Create annotation samplesheet
    echo -e "${annotation_csv}" > annotation.csv

    preprocess_rds.R "${sample}" "${rds_path}" meta.csv annotation.csv
    """
}