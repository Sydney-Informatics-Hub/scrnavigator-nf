// Load modules
include { ANNOTATE_CELL_CYCLE } from '../modules/annotate_cell_cycle'
include { ANNOTATE_DATABASE } from '../modules/annotate_db'

workflow ANNOTATE {
    take:
    integrated_rds
    species_params

    main:
    // Perform cell cycle annotation if required annotation data is provided
    valid_species_params = species_params.filter { p -> {
        p.species == 'human' || (
            p.s_genes != null &&
            p.g2m_genes != null
        )
    } }
    ANNOTATE_CELL_CYCLE(integrated_rds, valid_species_params)

    // Perform database-based annotation if required annotation data is provided
    db_annotation_in = ANNOTATE_CELL_CYCLE.out.annotated_rds
        .ifEmpty(integrated_rds)
    ANNOTATE_DATABASE(db_annotation_in, species_params)

    // Perform custom cell type annotation if required annotation data is provided
    custom_annotation_in = ANNOTATE_DATABASE.out.annotated_rds
        .ifEmpty(db_annotation_in)
    // ANNOTATE_CUSTOM(custom_annotation_in, species_params)  // TODO

    // Annotate clusters based on majority cell type
    // cluster_annotation_in = ANNOTATE_CUSTOM.out.annotated_rds
    //     .ifEmpty(custom_annotation_in)
    // ANNOTATE_CLUSTERS(cluster_annotation_in, species_params)  // TODO

    emit:
    rds = ANNOTATE_DATABASE.out.annotated_rds
    qc_results = ANNOTATE_DATABASE.out.qc_results
}