// Load modules
include { INTEGRATION } from '../modules/integrate'
include { CLUSTER_INTEGRATED } from '../modules/cluster_integrated'

workflow INTEGRATE {
    take:
    all_rds_to_integrate
    cohort_id
    integration_params

    main:
    // Perform integration
    INTEGRATION(all_rds_to_integrate, cohort_id)

    // Cluster integrated data
    cluster_input = INTEGRATION.out.integrated_rds
        .join(integration_params, by: 0)
    CLUSTER_INTEGRATED(cluster_input)

    emit:
    integrated_rds = CLUSTER_INTEGRATED.out.integrated_rds
    qc_results = CLUSTER_INTEGRATED.out.qc_results
}