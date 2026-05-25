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
    ens_db
    annotate_mt
    mt_gene_list
    vars_to_regress

    main:
    // Pre-process RDS files
    rds_files = samplesheet
        .map { row -> [ row.sample, row.rds_path, row.meta ] }
        .merge(species)
        .merge(ens_db)
        .merge(annotate_mt)
        .merge(mt_gene_list)
        .map { smp, rds, meta, spc, ens, mtann, mtlist -> {
            def ens_opt = ens == null ? [] : ens
            def mtlist_opt = mtlist == null ? [] : mtlist
            return [ smp, rds, meta, spc, ens_opt, mtann, mtlist_opt ]
        } }

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
        .map { row -> [ row.sample, row.params.res ]}
        .merge(all_resolutions)
        .merge(cluster_method)
    sct_in = FILTER.out.qc_rds
        .join(cluster_params, by: 0)
        .merge(vars_to_regress)
    SCTRANSFORM(sct_in)

    emit:
    rds = SCTRANSFORM.out.sct_rds
    filter_qc_results = FILTER.out.qc_results
    cluster_qc_results = SCTRANSFORM.out.qc_results
}