process PREPROCESS_RDS {
    publishDir "${params.outdir}/qc/${sample}/preprocess", mode: 'copy'

    input:
    tuple val(sample), path(rds_path), val(meta)
    val species_params

    output:
    tuple val(sample), path("${sample}.preprocessed.rds"), emit: preprocessed_rds

    script:
    def meta_csv = 'field,value\n' +
        meta.collect { field, value -> [ field, ( value == null ? '' : value ) ].join(',') }.join('\n')
    def species_csv = 'param,value\n' +
        species_params.collect { param, value -> [ param, ( value == null ? '' : value ) ].join(',') }.join('\n')
    """
    # Create parameter samplesheet
    echo -e "${meta_csv}" > meta.csv

    # Create species samplesheet
    echo -e "${species_csv}" > species.csv

    preprocess_rds.R "${sample}" "${rds_path}" meta.csv species.csv
    """
}