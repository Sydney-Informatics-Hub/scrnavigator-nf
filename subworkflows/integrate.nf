// Load modules
include { INTEGRATION } from '../modules/integrate'
include { CLUSTER_INTEGRATED } from '../modules/cluster_integrated'

workflow INTEGRATE {
    take:
    all_rds_to_integrate
    cohort_id
    all_resolutions
    cluster_method
    integrated_resolution

    main:
    // Perform integration
    INTEGRATION(all_rds_to_integrate, cohort_id)

    // Cluster integrated data
    cluster_input = INTEGRATION.out.integrated_rds
        .merge(all_resolutions)
        .merge(cluster_method)
        .merge(integrated_resolution)
    CLUSTER_INTEGRATED(cluster_input)

    emit:
    integrated_rds = CLUSTER_INTEGRATED.out.integrated_rds
    qc_results = INTEGRATION.out.qc_results
    cluster_qc_results = CLUSTER_INTEGRATED.out.qc_results
    gene_symbols = INTEGRATION.out.gene_symbols
    clustree_version = CLUSTER_INTEGRATED.out.version
}