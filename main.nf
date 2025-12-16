// Load modules
include { QUALITY_CONTROL } from './subworkflows/qc'
include { DETECT_DOUBLETS } from './modules/doublet'
include { INTEGRATE } from './subworkflows/integrate'
include { ANNOTATE } from './subworkflows/annotate'

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

    // Check annotation-related parameters are set
    assert params.species : "Error: Must provide species name."
    if (!params.no_mt) {
        assert ['human', 'mouse'].contains(params.species.toLowerCase()) ||
            !!params.mt_gene_list :
            "Error: If --species is neither 'human' nor 'mouse', --mt_gene_list must be provided."
    }
    assert ['human', 'mouse'].contains(params.species.toLowerCase()) ||
        !!params.ens_db_rds :
        "Error: If --species is neither 'human' nor 'mouse', --ens_db_rds must be provided."

    // Create channels from params
    // Required parameters and parameters with defaults
    all_resolutions                 = channel.value(params.resolutions)
    cluster_method                  = channel.value(params.cluster_method)
    integrated_resolution           = channel.value(params.integrated_resolution)
    cohort_id                       = channel.value(params.cohort_id)
    species                         = channel.value(params.species.toLowerCase())
    annotate_mt                     = channel.value(!params.no_mt)
    min_cells_for_annotation        = channel.value(params.min_cells_for_annotation as Integer)
    custom_annotation_mad_threshold = channel.value(params.custom_annotation_mad_threshold as Float)
    cell_type_proportion_threshold  = channel.value(params.cell_type_proportion_threshold as Float)

    // Optional parameters
    cluster_annotation         = !!params.cluster_annotation         ? channel.value(params.cluster_annotation)                                         : channel.value([null])
    mt_gene_list               = !!params.mt_gene_list               ? channel.fromPath(params.mt_gene_list, checkIfExists: true).first()               : channel.value([null])
    ens_db_rds                 = !!params.ens_db_rds                 ? channel.fromPath(params.ens_db_rds, checkIfExists: true).first()                 : channel.value([null])
    s_genes                    = !!params.s_genes                    ? channel.fromPath(params.s_genes, checkIfExists: true).first()                    : channel.value([null])
    g2m_genes                  = !!params.g2m_genes                  ? channel.fromPath(params.g2m_genes, checkIfExists: true).first()                  : channel.value([null])
    annotation_db              = !!params.annotation_db              ? channel.fromPath(params.annotation_db, checkIfExists: true).first()              : channel.value([null])
    custom_marker_genes        = !!params.custom_marker_genes        ? channel.fromPath(params.custom_marker_genes, checkIfExists: true).first()        : channel.value([null])
    manual_cluster_annotations = !!params.manual_cluster_annotations ? channel.fromPath(params.manual_cluster_annotations, checkIfExists: true).first() : channel.value([null])

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
            def cells_to_remove = row.cells_to_remove ? file(row.cells_to_remove, checkIfExists: true) : []

            def sample_params = [
                res:row.res,
                multiplet_rate:row.multiplet_rate,
                min_ncount:min_ncount,
                max_ncount:max_ncount,
                min_nfeature:min_nfeature,
                max_nfeature:max_nfeature,
                min_mt_pct:min_mt_pct,
                max_mt_pct:max_mt_pct
            ]
            def skip_keys = [ 'sample', 'rds', 'cells_to_remove' ] + sample_params.collect { key, _value -> key }
            def sample_meta = row.findAll { key, _value -> !skip_keys.contains(key) }
            return [
                sample:sample,
                rds_path:rds_path,
                params:sample_params,
                meta:sample_meta,
                cells_to_remove:cells_to_remove
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
    QUALITY_CONTROL(
        samplesheet,
        all_resolutions,
        cluster_method,
        species,
        ens_db_rds,
        annotate_mt,
        mt_gene_list
    )

    // If only running QC, stop here, otherwise continue on
    if (!params.qc_only) {
        // Doublet detection
        doublet_params = samplesheet
            .map { row -> [ row.sample, row.multiplet_rate, row.res ]}
            .merge(all_resolutions)
            .merge(cluster_method)
        doublet_in = QUALITY_CONTROL.out.rds
            .join(doublet_params, by: 0)
        DETECT_DOUBLETS(doublet_in)

        // Integration
        all_rds_to_integrate = DETECT_DOUBLETS.out.doublets_removed_rds
            .map { _sample, rds -> rds }
            .collect()
        INTEGRATE(
            all_rds_to_integrate,
            cohort_id,
            all_resolutions,
            cluster_method,
            integrated_resolution
        )
    }

    if (!params.qc_only && !params.no_analysis) {
        annotation_rds_input = INTEGRATE.out.integrated_rds
        ANNOTATE(
            annotation_rds_input,
            species,
            s_genes,
            g2m_genes,
            ens_db_rds,
            min_cells_for_annotation,
            annotation_db,
            custom_marker_genes,
            custom_annotation_mad_threshold,
            cluster_annotation,
            cell_type_proportion_threshold,
            manual_cluster_annotations
        )
        // TODO:
        // PSEUDO()
        // DE()
        // FEA
    }

}