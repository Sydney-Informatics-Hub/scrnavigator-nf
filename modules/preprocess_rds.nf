process PREPROCESS_RDS {
    publishDir "${params.outdir}/qc/${sample}/preprocess", mode: 'copy'
    ext input_size: { new InputFileSizes(rds_path) + new InputFileSizes(ens_db_rds) }
    memory { task.ext.input_size.getSizeMB() * 2 + 1.GB }

    input:
    tuple val(sample), path(rds_path), val(meta), val(species), path(ens_db_rds), val(annotate_mt), path(mt_gene_list)

    output:
    tuple val(sample), path("${sample}.preprocessed.rds"), emit: preprocessed_rds

    script:
    def meta_csv = 'field,value\n' +
        meta.collect { field, value -> [ field, ( value == null ? '' : value ) ].join(',') }.join('\n')
    def annotation_csv = 'param,value\n' +
        "species,${species}\n" +
        "ens_db_rds,${ens_db_rds}\n" +
        "annotate_mt,${annotate_mt}\n" +
        "mt_gene_list,${mt_gene_list}\n"
    """
    # Create parameter samplesheet
    echo -e "${meta_csv}" > meta.csv

    # Create annotation samplesheet
    echo -e "${annotation_csv}" > annotation.csv

    preprocess_rds.R "${sample}" "${rds_path}" meta.csv annotation.csv
    """
}