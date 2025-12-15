// Load modules
include { QUALITY_CONTROL } from './subworkflows/qc'
include { DETECT_DOUBLETS } from './modules/doublet'
include { INTEGRATE } from './subworkflows/integrate'

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

    // Check species-related parameters are set
    assert params.species : "Error: Must provide species name."
    if (!params.no_mt) {
        assert ['human', 'mouse'].contains(params.species.toLowerCase()) ||
            !!params.mt_gene_list :
            "Error: If --species is neither 'human' nor 'mouse', --mt_gene_list must be provided."
    }
    assert ['human', 'mouse'].contains(params.species.toLowerCase()) ||
            !!params.ens_db_rds :
            "Error: If --species is neither 'human' nor 'mouse', --ens_db_rds must be provided."
    species_params = channel.value([
        species:params.species.toLowerCase(),
        annotate_mt:!params.no_mt,
        mt_gene_list:(!!params.mt_gene_list ? file(params.mt_gene_list, checkIfExists: true) : null),
        ens_db_rds:(!!params.ens_db_rds ? file(params.ens_db_rds, checkIfExists: true) : null)
    ])

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
            def cells_to_remove = row.cells_to_remove ? file(row.cells_to_remove, checkIfExists: true) : null

            def sample_params = [
                res:row.res,
                cells_to_remove:cells_to_remove,
                multiplet_rate:row.multiplet_rate,
                min_ncount:min_ncount,
                max_ncount:max_ncount,
                min_nfeature:min_nfeature,
                max_nfeature:max_nfeature,
                min_mt_pct:min_mt_pct,
                max_mt_pct:max_mt_pct
            ]
            def skip_keys = [ 'sample', 'rds' ] + sample_params.collect { key, _value -> key }
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

    // Run initial quality control
    all_resolutions = channel.value(params.resolutions)
    cluster_method = channel.value(params.cluster_method)

    QUALITY_CONTROL(samplesheet, all_resolutions, cluster_method, species_params)

    // If only running QC, stop here, otherwise continue on
    if (!params.qc_only) {
        // Doublet detection
        DETECT_DOUBLETS(QUALITY_CONTROL.out.rds, all_resolutions, cluster_method)

        // Integration
        all_rds_to_integrate = DETECT_DOUBLETS.out.doublets_removed_rds
            .map { _sample, rds, _meta -> rds }
            .collect()
        integration_params = channel.of([
            params.cohort_id,
            [
                resolutions:params.resolutions,
                integrated_resolution:params.integrated_resolution,
                cluster_method:params.cluster_method
            ]
        ])
        cohort_id = channel.value(params.cohort_id)
        INTEGRATE(all_rds_to_integrate, cohort_id, integration_params)
    }

    if (!params.qc_only && !params.no_analysis) {
        // TODO:
        // ANNOTATE()
        // PSEUDO()
        // DE()
        // FEA
    }

}