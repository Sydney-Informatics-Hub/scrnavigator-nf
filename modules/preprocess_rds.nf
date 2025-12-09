process PREPROCESS_RDS {
    input:
    tuple val(sample), path(rds_path), val(ensembl)

    output:
    tuple val(sample), path("${sample}.preprocessed.rds"), emit: preprocessed_rds

    script:
    def is_ensembl = ensembl ? "TRUE" : "FALSE"
    """
    preprocess_rds.R "${sample}" "${rds_path}" ${is_ensembl}
    """
}