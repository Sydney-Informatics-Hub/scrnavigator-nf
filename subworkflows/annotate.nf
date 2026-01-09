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
    valid_cc_annotation_params = species
        .merge(s_genes)
        .merge(g2m_genes)
        .merge(ens_db_rds)
        .filter { spc, sg, gg, _ens -> {
            spc == 'human' || (
                sg != null &&
                gg != null
            )
        } }
        .map { spc, sg, gg, ens -> {
            def sg_opt = sg == null ? [] : sg
            def gg_opt = gg == null ? [] : gg
            def ens_opt = ens == null ? [] : ens
            return [ spc, sg_opt, gg_opt, ens_opt ]
         } }
    cc_annotation_in = integrated_rds
        .merge(valid_cc_annotation_params)

    ANNOTATE_CELL_CYCLE(cc_annotation_in)

    cc_annotation_out = ANNOTATE_CELL_CYCLE.out.annotated_rds
        .join(integrated_rds, by: 0, remainder: true)
        .map { cohort, cc_rds, int_rds -> {
            def rds = cc_rds != null ? cc_rds : int_rds
            return [ cohort, rds ]
        } }

    // Perform database-based annotation if required annotation data is provided
    valid_db_annotation_params = species
        .merge(min_cells_for_annotation)
        .merge(annotation_db)
        .filter { spc, _mc, anndb -> {
            spc == 'human' || anndb != null
        } }
        .map { spc, mc, anndb -> {
            def anndb_opt = anndb == null ? [] : anndb
            return [ spc, mc, anndb_opt ]
        } }
    db_annotation_in = cc_annotation_out
        .merge(valid_db_annotation_params)

    ANNOTATE_DATABASE(db_annotation_in)

    db_annotation_out = ANNOTATE_DATABASE.out.annotated_rds
        .join(cc_annotation_out, by: 0, remainder: true)
        .map { cohort, db_rds, cc_rds -> {
            def rds = db_rds != null ? db_rds : cc_rds
            return [ cohort, rds ]
        } }

    // Perform custom cell type annotation if required annotation data is provided
    valid_custom_annotation_params = species
        .merge(ens_db_rds)
        .merge(custom_marker_genes)
        .merge(custom_annotation_mad_threshold)
        .filter { _spc, _ens, cus, _mad -> {
            cus != null
        } }
        .map { spc, ens, cus, mad -> {
            def ens_opt = ens == null ? [] : ens
            def cus_opt = cus == null ? [] : cus
            return [ spc, ens_opt, cus_opt, mad ]
        } }
    custom_annotation_in = db_annotation_out
        .merge(valid_custom_annotation_params)

    ANNOTATE_CUSTOM(custom_annotation_in)

    custom_annotation_out = ANNOTATE_CUSTOM.out.annotated_rds
        .join(ANNOTATE_DATABASE.out.annotated_rds, by: 0, remainder: true)
        .join(ANNOTATE_CELL_CYCLE.out.annotated_rds, by: 0, remainder: true)
        .map { cohort, cus_rds, db_rds, cc_rds -> {
            def rds = cus_rds != null ? cus_rds : ( db_rds != null ? db_rds : cc_rds )
            return [ cohort, rds ]
        } }
        .filter { _cohort, rds -> rds != null }

    // Annotate clusters based on majority cell type
    // Runs if any annotation has been performed
    valid_cluster_annotation_params = species
        .merge(cluster_annotation)
        .merge(cell_type_proportion_threshold)
        .merge(manual_cluster_annotations)
        .map { spc, clu, ctprop, annfile -> {
            def annfile_opt = annfile == null ? [] : annfile
            return [ spc, clu, ctprop, annfile_opt ]
        } }
    cluster_annotation_in = custom_annotation_out
        .merge(valid_cluster_annotation_params)

    ANNOTATE_CLUSTERS(cluster_annotation_in)

    // Return cluster-annotated RDS if any annotation was performed,
    // otherwise return original integrated RDS
    annotated_rds = ANNOTATE_CLUSTERS.out.annotated_rds
        .join(integrated_rds, by: 0, remainder: true)
        .map { cohort, clu_rds, int_rds -> {
            def rds = clu_rds != null ? clu_rds : int_rds
            return [ cohort, rds ]
        } }

    emit:
    rds = annotated_rds
    cell_cycle_annotation_qc_results = ANNOTATE_CELL_CYCLE.out.qc_results
    db_annotation_qc_results = ANNOTATE_DATABASE.out.qc_results
    custom_annotation_qc_results = ANNOTATE_CUSTOM.out.qc_results
    cluster_annotation_qc_results = ANNOTATE_CLUSTERS.out.qc_results
}