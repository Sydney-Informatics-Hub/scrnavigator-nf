// Load modules
include { PREPROCESS_RDS } from '../modules/preprocess_rds'
include { FILTER } from '../modules/filter'
include { SCTRANSFORM } from '../modules/sct'

workflow QUALITY_CONTROL {
    take:
    samplesheet
    all_resolutions
    cluster_method
    species
    ens_db_rds
    annotate_mt
    mt_gene_list

    main:
    // Pre-process RDS files
    rds_files = samplesheet
        .map { row -> [ row.sample, row.rds_path, row.meta ] }
        .merge(species)
        .merge(ens_db_rds)
        .merge(annotate_mt)
        .merge(mt_gene_list)

    // Conduct initial QC
    PREPROCESS_RDS(rds_files)

    // Perform filtering
    sample_parameters = samplesheet
        .map { row -> [ row.sample, row.params, row.cells_to_remove ] }
    filter_in = PREPROCESS_RDS.out.preprocessed_rds
        .join(sample_parameters, by: 0)

    FILTER(filter_in)

    // Run SCTransform
    cluster_params = samplesheet
        .map { row -> [ row.sample, row.res ]}
        .merge(all_resolutions)
        .merge(cluster_method)
    sct_in = FILTER.out.qc_rds
        .join(cluster_params, by: 0)
    SCTRANSFORM(sct_in)

    emit:
    rds = SCTRANSFORM.out.sct_rds
    filter_qc_results = FILTER.out.qc_results
    cluster_qc_results = SCTRANSFORM.out.qc_results
}