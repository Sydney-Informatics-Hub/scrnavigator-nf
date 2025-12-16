// Load modules
include { ANNOTATE_CELL_CYCLE } from '../modules/annotate_cell_cycle'
include { ANNOTATE_DATABASE } from '../modules/annotate_db'
include { ANNOTATE_CUSTOM } from '../modules/annotate_custom'
include { ANNOTATE_CLUSTERS } from '../modules/annotate_clusters'

workflow ANNOTATE {
    take:
    integrated_rds
    species
    s_genes
    g2m_genes
    ens_db_rds
    min_cells_for_annotation
    annotation_db
    custom_marker_genes
    custom_annotation_mad_threshold
    cluster_annotation
    cell_type_proportion_threshold
    manual_cluster_annotations

    main:
    // Perform cell cycle annotation if required annotation data is provided
    cc_annotation_params = species
        .merge(s_genes)
        .merge(g2m_genes)
        .merge(ens_db_rds)
    valid_cc_annotation_params = cc_annotation_params
        .filter { spc, sg, gg, _ens -> {
            spc == 'human' || (
                sg != null &&
                gg != null
            )
        } }
        .map { spc, sg, gg, ens -> [
            spc,
            (sg != null ? sg : []),
            (gg != null ? gg : []),
            (ens != null ? ens : [])
        ] }
    ANNOTATE_CELL_CYCLE(integrated_rds, valid_cc_annotation_params)

    // Perform database-based annotation if required annotation data is provided
    valid_db_annotation_params = annotation_params.filter { p -> {
        p.species == 'human' || (
            p.annotation_db != null
        )
    } }
    db_annotation_in = ANNOTATE_CELL_CYCLE.out.annotated_rds
        .ifEmpty(integrated_rds)
    ANNOTATE_DATABASE(db_annotation_in, valid_db_annotation_params)

    // Perform custom cell type annotation if required annotation data is provided
    valid_custom_annotation_params = annotation_params.filter { p -> {
        p.custom_marker_genes != null
    } }
    custom_annotation_in = ANNOTATE_DATABASE.out.annotated_rds
        .ifEmpty(db_annotation_in)
    ANNOTATE_CUSTOM(custom_annotation_in, valid_custom_annotation_params)

    // Annotate clusters based on majority cell type
    cluster_annotation_in = ANNOTATE_CUSTOM.out.annotated_rds
        .ifEmpty(ANNOTATE_DATABASE.out.annotated_rds)
        .ifEmpty(ANNOTATE_CELL_CYCLE.out.annotated_rds)
    ANNOTATE_CLUSTERS(cluster_annotation_in, annotation_params)

    emit:
    rds = ANNOTATE_CLUSTERS.out.annotated_rds
    cell_cycle_annotation_qc_results = ANNOTATE_CELL_CYCLE.out.qc_results
    db_annotation_qc_results = ANNOTATE_DATABASE.out.qc_results
    custom_annotation_qc_results = ANNOTATE_CUSTOM.out.qc_results
    cluster_annotation_qc_results = ANNOTATE_CLUSTERS.out.qc_results
}