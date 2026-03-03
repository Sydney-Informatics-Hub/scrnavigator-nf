process ANNOTATE_CELL_CYCLE {
    publishDir "${params.outdir}/annotation/cell_cycle", mode: 'copy'
    // ext input_size: { new InputFileSizes(rds_path) + new InputFileSizes(ens_db_rds) }
    ext input_size: { rds_path.size() + ( ens_db_rds ? ens_db_rds.size() : 0 ) }
    memory { 8.GB + task.ext.input_size.B * 10 }

    input:
    tuple val(cohort_name), path(rds_path), val(species), path(s2_genes), path(g2m_genes), path(ens_db_rds)

    output:
    tuple val(cohort_name), path("${cohort_name}.annotated.cell_cycle.rds"), emit: annotated_rds
    tuple val(cohort_name), path("available_annotations.txt"), emit: available_annotations
    tuple val(cohort_name), path("qc_results"), emit: qc_results

    script:
    def annotation_csv = 'param,value\n' +
        "species,${species}\n" +
        "s2_genes,${s2_genes}\n" +
        "g2m_genes,${g2m_genes}\n" +
        "ens_db_rds,${ens_db_rds}\n"
    """
    # Create annotation samplesheet
    echo -e "${annotation_csv}" > annotation.csv

    annotate_cell_cycle.R "${cohort_name}" "${rds_path}" annotation.csv
    """
}