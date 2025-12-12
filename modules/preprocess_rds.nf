process PREPROCESS_RDS {
    input:
    tuple val(sample), path(rds_path), val(meta)

    output:
    tuple val(sample), path("${sample}.preprocessed.rds"), emit: preprocessed_rds

    script:
    def meta_csv = 'field,value\n' +
        meta.collect { field, value -> [ field, ( value ? value : '' ) ].join(',') }.join('\n')
    """
    # Create parameter samplesheet
    echo -e "${meta_csv}" > meta.csv

    preprocess_rds.R "${sample}" "${rds_path}" meta.csv
    """
}