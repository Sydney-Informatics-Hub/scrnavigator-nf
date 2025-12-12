// Load modules
include { PREPROCESS_RDS } from '../modules/preprocess_rds'
include { INIT_QC } from '../modules/init_qc'

workflow QUALITY_CONTROL {
    take:
    samplesheet
    all_resolutions
    cluster_method

    main:
    // Pre-process RDS files
    rds_files = samplesheet
        .map { row -> [ row.sample, row.rds_path, row.meta ] }

    PREPROCESS_RDS(rds_files)

    // Conduct initial QC
    sample_parameters = samplesheet
        .map { row -> [ row.sample, row.params ] }
    init_qc_in = PREPROCESS_RDS.out.preprocessed_rds
        .join(sample_parameters, by: 0)

    INIT_QC(init_qc_in, all_resolutions, cluster_method)

    emit:
    qc_data = INIT_QC.out.qc_rds
    qc_results = INIT_QC.out.qc_results
}