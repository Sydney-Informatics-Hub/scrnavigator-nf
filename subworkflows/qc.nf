// Load modules
include { PREPROCESS_RDS } from '../modules/preprocess_rds'
include { FILTER } from '../modules/filter'
include { SCTRANSFORM } from '../modules/sct'

workflow QUALITY_CONTROL {
    take:
    samplesheet
    all_resolutions
    cluster_method
    annotation_params

    main:
    // Pre-process RDS files
    rds_files = samplesheet
        .map { row -> [ row.sample, row.rds_path, row.meta ] }

    PREPROCESS_RDS(rds_files, annotation_params)

    // Conduct initial QC
    sample_parameters = samplesheet
        .map { row -> [ row.sample, row.params ] }
    filter_in = PREPROCESS_RDS.out.preprocessed_rds
        .join(sample_parameters, by: 0)

    FILTER(filter_in)

    SCTRANSFORM(FILTER.out.qc_rds, all_resolutions, cluster_method)

    emit:
    rds = SCTRANSFORM.out.sct_rds
    filter_qc_results = FILTER.out.qc_results
    cluster_qc_results = SCTRANSFORM.out.qc_results
}