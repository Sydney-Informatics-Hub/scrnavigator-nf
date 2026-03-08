process INTEGRATION {
    publishDir "${params.outdir}/integration", mode: 'copy'
    // ext input_size: { new InputFileSizes(rds_paths) }
    ext input_size: { rds_paths instanceof Collection ? rds_paths.inject(0) { sum, f -> sum + f.size() } : rds_paths.size() }
    memory { 16.GB + task.ext.input_size.B * 10 }
    container "sydneyinformaticshub/scrnavigator-nf-cluster"

    input:
    path rds_paths
    val cohort_name

    output:
    tuple val(cohort_name), path("${cohort_name}.integrated.rds"), emit: integrated_rds
    tuple val(cohort_name), path("qc_results"), emit: qc_results
    tuple val(cohort_name), path("gene_symbols.Rds"), emit: gene_symbols

    script:
    def all_rds_paths = rds_paths.collect { f -> "'${f}'" }.join(' ')
    """
    integrate.R "${cohort_name}" ${all_rds_paths}
    """
}