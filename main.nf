// Load modules
include { PREPROCESS_RDS } from './modules/preprocess_rds.nf'
include { INIT_QC } from './modules/init_qc.nf'

// Required pipeline parameters
params.help = false
params.input = null

def helpMessage() {
    log.info"""
    Usage: nextflow run main.nf --input /path/to/samplesheet.csv 

    Required parameters:

    --input               Spectify full path and name of samplesheet csv

    Optional parameters:

    --outdir              Specify path to output directory.

    """.stripIndent()
}

workflow {
    // Print help message if requested or if no input is provided
    if (params.help || !params.input) {
        helpMessage()
        System.exit(1)
    }

    // Read in samplesheet
    samplesheet = channel.fromPath(params.input)
        .splitCsv( header: true )
        .map { row -> {
            def sample = row.sample
            assert !!sample : "Error: Must provide a sample name in the samplesheet."
            def rds_path = file(row.rds, checkIfExists: true)
            assert !!rds_path : "Error: Must provide a valid path to a Seurat RDS file."
            def min_ncount = row.min_ncount ? row.min_ncount.toInteger() : 0
            def max_ncount = row.max_ncount ? row.max_ncount.toInteger() : null
            def min_nfeature = row.min_nfeature ? row.min_nfeature.toInteger() : 0
            def max_nfeature = row.max_nfeature ? row.max_nfeature.toInteger() : null
            def min_mt_pct = row.min_mt_pct ? row.min_mt_pct.toInteger() : 0
            def max_mt_pct = row.max_mt_pct ? row.max_mt_pct.toInteger() : 100

            def sample_params = [
                res:row.res,
                clusters_to_remove:row.clusters_to_remove,
                multiplet_rate:row.multiplet_rate,
                min_ncount:min_ncount,
                max_ncount:max_ncount,
                min_nfeature:min_nfeature,
                max_nfeature:max_nfeature,
                min_mt_pct:min_mt_pct,
                max_mt_pct:max_mt_pct
            ]
            def skip_keys = [ 'sample', 'rds_path' ] + sample_params.collect { key, _value -> key }
            def sample_meta = row.findAll { key, _value -> !skip_keys.contains(key) }
            return [
                sample:sample,
                rds_path:rds_path,
                params:sample_params,
                meta:sample_meta
            ]
        }}

    // Validate samplesheet: only one entry per sample
    samplesheet
        .map { row -> [ row.sample, row.rds_path ] }
        .groupTuple()
        .map { sample, rds_paths -> {
            assert rds_paths.size() != 0 : "Error: Missing RDS file for sample '${sample}'."
            assert rds_paths.size() == 1 : "Error: Duplicate entries found for sample '${sample}'."
            return [ sample, rds_paths ]
        }}

    // Pre-process RDS files
    rds_files = samplesheet
        .map { row -> [ row.sample, row.rds_path, row.meta ] }

    PREPROCESS_RDS(rds_files)

    // Conduct initial QC
    sample_parameters = samplesheet
        .map { sample_id, _rds_path, sample_params, _sample_meta -> [ sample_id, sample_params ] }
    init_qc_in = PREPROCESS_RDS.out.preprocessed_rds
        .join(sample_parameters, by: 0)
    all_resolutions = channel.value(params.resolutions)
    cluster_method = channel.value(params.cluster_method)

    INIT_QC(init_qc_in, all_resolutions, cluster_method)


}