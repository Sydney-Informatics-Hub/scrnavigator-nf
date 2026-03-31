// Load modules
include { QUALITY_CONTROL } from './subworkflows/qc'
include { DETECT_DOUBLETS } from './modules/doublet'
include { INTEGRATE } from './subworkflows/integrate'
include { ANNOTATE } from './subworkflows/annotate'
include { PSEUDOBULK } from './modules/pseudobulk'
include { DIFFERENTIAL_EXPRESSION } from './modules/de'
include { ANALYSE_DE } from './modules/analyse_de'
include { ORA } from './modules/ora'
include { GSEA } from './modules/gsea'
include { ANALYSE_ORA } from './modules/analyse_ora'
include { ANALYSE_GSEA } from './modules/analyse_gsea'
include { REPORT } from './modules/report'

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
    p_value_threshold               = channel.value(params.p_value_threshold as Float)

    // Optional parameters
    cluster_annotation         = !!params.cluster_annotation         ? channel.value(params.cluster_annotation)                                         : channel.value([null])
    mt_gene_list               = !!params.mt_gene_list               ? channel.fromPath(params.mt_gene_list, checkIfExists: true).first()               : channel.value([null])
    ens_db_rds                 = !!params.ens_db_rds                 ? channel.fromPath(params.ens_db_rds, checkIfExists: true).first()                 : channel.value([null])
    s_genes                    = !!params.s_genes                    ? channel.fromPath(params.s_genes, checkIfExists: true).first()                    : channel.value([null])
    g2m_genes                  = !!params.g2m_genes                  ? channel.fromPath(params.g2m_genes, checkIfExists: true).first()                  : channel.value([null])
    annotation_db              = !!params.annotation_db              ? channel.fromPath(params.annotation_db, checkIfExists: true).first()              : channel.value([null])
    custom_marker_genes        = !!params.custom_marker_genes        ? channel.fromPath(params.custom_marker_genes, checkIfExists: true).first()        : channel.value([null])
    manual_cluster_annotations = !!params.manual_cluster_annotations ? channel.fromPath(params.manual_cluster_annotations, checkIfExists: true).first() : channel.value([null])
    pseudo_groups              = !!params.pseudo_groups              ? channel.value(params.pseudo_groups)                                              : channel.value([null])
    fc_threshold               = !!params.fc_threshold               ? channel.value(params.fc_threshold as Float)                                      : channel.value([null])
    ora_db                     = !!params.ora_db                     ? channel.fromPath(params.ora_db, checkIfExists: true).first()                     : channel.value([])
    gsea_db                    = !!params.gsea_db                    ? channel.fromPath(params.gsea_db, checkIfExists: true).first()                    : channel.value([])

    // Read in samplesheet
    samplesheet_csv = channel.fromPath(params.input, checkIfExists: true)
    samplesheet = samplesheet_csv
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
            .map { row -> [ row.sample, row.params.multiplet_rate, row.params.res ]}
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

        detect_doublets_qc = DETECT_DOUBLETS.out.qc_results
        integrate_qc = INTEGRATE.out.qc_results
        integrate_cluster_qc = INTEGRATE.out.cluster_qc_results
    } else {
        detect_doublets_qc = channel.empty()
        integrate_qc = channel.empty()
        integrate_cluster_qc = channel.empty()
    }

    if (!params.qc_only && !params.no_analysis) {
        // Annotation
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

        // Pseudobulking
        pseudo_in = ANNOTATE.out.rds
            .merge(pseudo_groups)
            .filter { _id, _rds, grps -> {
                grps != null
            } }
        PSEUDOBULK(pseudo_in)

        // Differential expression
        comparisons = !params.comparisons ? channel.empty() : (
            channel.fromPath(params.comparisons)
                .splitCsv( header: true )
                .map { row -> {
                    assert row.ref != null && row.ref != ''
                    assert row.test != null && row.test != ''
                    [ row.ref, row.test ]
                } }
        )
        de_in = PSEUDOBULK.out.pseudobulked_rds
            .combine(comparisons)
        DIFFERENTIAL_EXPRESSION(de_in)

        // Merge and analyse DE results
        all_de_results = DIFFERENTIAL_EXPRESSION.out.de_rds
            .groupTuple()
            .join(INTEGRATE.out.gene_symbols)
            .merge(p_value_threshold)
            .merge(fc_threshold)
        ANALYSE_DE(all_de_results)

        // Functional enrichment analysis
        fea_in = ANALYSE_DE.out.de_rds
            .combine(comparisons)
            .merge(species)
        ORA(fea_in, ora_db)
        GSEA(fea_in, gsea_db)

        // Merge and analyse FEA results
        all_ora_results = ORA.out.ora_rds
            .groupTuple()
        ANALYSE_ORA(all_ora_results)
        all_gsea_results = GSEA.out.gsea_rds
            .groupTuple()
        ANALYSE_GSEA(all_gsea_results)

        annotate_cc_qc = ANNOTATE.out.cell_cycle_annotation_qc_results
        annotate_db_qc = ANNOTATE.out.db_annotation_qc_results
        annotate_custom_qc = ANNOTATE.out.custom_annotation_qc_results
        annotate_cluster_qc = ANNOTATE.out.cluster_annotation_qc_results
        pseudo_comp = PSEUDOBULK.out.comparison_groups
        de_csv = ANALYSE_DE.out.de_csv
        de_results = ANALYSE_DE.out.de_results
        gsea_csv = ANALYSE_GSEA.out.gsea_csv
        gsea_reduced_csv = ANALYSE_GSEA.out.gsea_reduced_csv
        gsea_plots = ANALYSE_GSEA.out.gsea_plots
        ora_csv = ANALYSE_ORA.out.ora_csv
        ora_reduced_csv = ANALYSE_ORA.out.ora_reduced_csv
        ora_plots = ANALYSE_ORA.out.ora_plots
        available_annotations = ANNOTATE.out.available_annotations
    } else {
        annotate_cc_qc = channel.empty()
        annotate_db_qc = channel.empty()
        annotate_custom_qc = channel.empty()
        annotate_cluster_qc = channel.empty()
        pseudo_comp = channel.empty()
        de_csv = channel.empty()
        de_results = channel.empty()
        gsea_csv = channel.empty()
        gsea_reduced_csv = channel.empty()
        gsea_plots = channel.empty()
        ora_csv = channel.empty()
        ora_reduced_csv = channel.empty()
        ora_plots = channel.empty()
        available_annotations = channel.empty()
    }

    // Summary report
    // Collect all results
    // QC - one or more 'qc_results' directories per stage
    // Collect each set of directories into a single tuple
    // If the stage wasn't run, create an empty nested tuple (for merging purposes)
    qc_filter_results = QUALITY_CONTROL.out.filter_qc_results
        .map { _smp, qc_dir -> qc_dir }
        .collect()
        .map { dirs -> [ dirs ] }
        .ifEmpty([[]])
    qc_cluster_results = QUALITY_CONTROL.out.cluster_qc_results
        .map { _smp, qc_dir -> qc_dir }
        .collect()
        .map { dirs -> [ dirs ] }
        .ifEmpty([[]])
    qc_doublet_results = detect_doublets_qc
        .map { _smp, qc_dir -> qc_dir }
        .collect()
        .map { dirs -> [ dirs ] }
        .ifEmpty([[]])
    // Integration
    // From now on, only one output directory per stage
    integration_qc_results = integrate_qc
        .map { _cohort, qc_dir -> qc_dir }
        .ifEmpty([[]])
    integration_cluster_results = integrate_cluster_qc
        .map { _cohort, qc_dir -> qc_dir }
        .ifEmpty([[]])
    // Annotation
    annotation_cc_results = annotate_cc_qc
        .map { _cohort, qc_dir -> qc_dir }
        .ifEmpty([[]])
    annotation_db_results = annotate_db_qc
        .map { _cohort, qc_dir -> qc_dir }
        .ifEmpty([[]])
    annotation_custom_results = annotate_custom_qc
        .map { _cohort, qc_dir -> qc_dir }
        .ifEmpty([[]])
    annotation_clusters_results = annotate_cluster_qc
        .map { _cohort, qc_dir -> qc_dir }
        .ifEmpty([[]])
    // Analysis
    analysis_pseudo_comparison_groups = pseudo_comp
        .map { _cohort, cmp_file -> cmp_file }
        .ifEmpty([[]])
    analysis_de_results = de_csv
        .mix(de_results)
        .map { _cohort, out_file -> out_file }
        .collect()
        .map { res -> [ res ] }
        .ifEmpty([[]])
    analysis_gsea_results = gsea_csv
        .mix(gsea_reduced_csv)
        .mix(gsea_plots)
        .map { _cohort, out_file -> out_file }
        .collect()
        .map { res -> [ res ] }
        .ifEmpty([[]])
    analysis_ora_results = ora_csv
        .mix(ora_reduced_csv)
        .mix(ora_plots)
        .map { _cohort, out_file -> out_file }
        .collect()
        .map { res -> [ res ] }
        .ifEmpty([[]])
    // Available annotation files
    available_annotation_files = available_annotations
        .map { _cohort, ann_file -> ann_file }
        .collect()
        .map { anns -> [ anns ] }
        .ifEmpty([[]])

    // Gather additional metadata
    report_metadata = samplesheet
        .map { row -> row.meta.keySet() }
        .flatten()
        .unique()
        .collect()
        .map { meta_fields -> [ meta_fields ] }
        .merge(integrated_resolution)
        .merge(pseudo_groups)
        .map { meta, int_res, pseudo -> [ meta_fields: meta, integration_resolution: int_res, pseudo_groups: pseudo ] }

    // Gather all report templates
    index_template = channel.fromPath("${projectDir}/assets/index.qmd", checkIfExists: true)
    qc_filter_template = channel.fromPath("${projectDir}/assets/qc_filter.qmd", checkIfExists: true)
        .merge(qc_filter_results)
        .filter { _t, r -> !!r }
        .map { t, _r -> t }
    qc_cluster_template = channel.fromPath("${projectDir}/assets/qc_cluster.qmd", checkIfExists: true)
        .merge(qc_cluster_results)
        .filter { _t, r -> !!r }
        .map { t, _r -> t }
    qc_doublet_template = channel.fromPath("${projectDir}/assets/qc_doublets.qmd", checkIfExists: true)
        .merge(qc_doublet_results)
        .filter { _t, r -> !!r }
        .map { t, _r -> t }
    integration_qc_template = channel.fromPath("${projectDir}/assets/integration_qc.qmd", checkIfExists: true)
        .merge(integration_qc_results)
        .filter { _t, r -> !!r }
        .map { t, _r -> t }
    integration_cluster_template = channel.fromPath("${projectDir}/assets/integration_cluster.qmd", checkIfExists: true)
        .merge(integration_cluster_results)
        .filter { _t, r -> !!r }
        .map { t, _r -> t }
    annotation_results = annotation_cc_results
        .mix(annotation_db_results)
        .mix(annotation_custom_results)
        .mix(annotation_clusters_results)
        .flatten()
        .collect()
        .map { res -> [ res ] }
    annotation_template = channel.fromPath("${projectDir}/assets/annotation.qmd", checkIfExists: true)
        .merge(annotation_results)
        .filter { _t, r -> !!r }
        .map { t, _r -> t }
    analysis_pseudo_template = channel.fromPath("${projectDir}/assets/analysis_pseudo.qmd", checkIfExists: true)
        .merge(analysis_pseudo_comparison_groups)
        .filter { _t, r -> !!r }
        .map { t, _r -> t }
    analysis_de_template = channel.fromPath("${projectDir}/assets/analysis_de.qmd", checkIfExists: true)
        .merge(analysis_de_results)
        .filter { _t, r -> !!r }
        .map { t, _r -> t }
    analysis_gsea_template = channel.fromPath("${projectDir}/assets/analysis_gsea.qmd", checkIfExists: true)
        .merge(analysis_gsea_results)
        .filter { _t, r -> !!r }
        .map { t, _r -> t }
    analysis_ora_template = channel.fromPath("${projectDir}/assets/analysis_ora.qmd", checkIfExists: true)
        .merge(analysis_ora_results)
        .filter { _t, r -> !!r }
        .map { t, _r -> t }
        .map { t, _r -> t }
    report_templates = index_template
        .mix(qc_filter_template)
        .mix(qc_cluster_template)
        .mix(qc_doublet_template)
        .mix(integration_qc_template)
        .mix(integration_cluster_template)
        .mix(annotation_template)
        .mix(analysis_pseudo_template)
        .mix(analysis_de_template)
        .mix(analysis_gsea_template)
        .mix(analysis_ora_template)
        .collect()
        .map { tmp -> [ tmp ] }
    report_style = channel.fromPath("${projectDir}/assets/styles.scss")

    // Gather CSV inputs
    custom_marker_genes_csv = !!params.custom_marker_genes ? channel.fromPath(params.custom_marker_genes, checkIfExists: true) : channel.value([])
    comparisons_csv = !!params.comparisons ? channel.fromPath(params.comparisons, checkIfExists: true) : channel.value([])

    // Run the report module
    report_input = cohort_id
        .merge(report_metadata)
        .merge(report_templates)
        .merge(qc_filter_results)
        .merge(qc_cluster_results)
        .merge(qc_doublet_results)
        .merge(integration_qc_results)
        .merge(integration_cluster_results)
        .merge(annotation_cc_results)
        .merge(annotation_db_results)
        .merge(annotation_custom_results)
        .merge(annotation_clusters_results)
        .merge(analysis_pseudo_comparison_groups)
        .merge(analysis_de_results)
        .merge(analysis_gsea_results)
        .merge(analysis_ora_results)
        .merge(available_annotation_files)
        .merge(report_style)
    REPORT(report_input, samplesheet_csv, custom_marker_genes_csv, comparisons_csv)
}