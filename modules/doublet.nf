process DETECT_DOUBLETS {
    publishDir "${params.outdir}/qc/${sample}/doublets", mode: 'copy'
    // ext input_size: { new InputFileSizes(rds_path) }
    ext input_size: { rds_path.size() }
    memory { 16.GB + task.ext.input_size.B * 10 }
    container "sydneyinformaticshub/scrnavigator-nf-doublet"

    input:
    tuple val(sample), path(rds_path, stageAs: "input/*"), val(multiplet_rate), val(default_res), val(all_resolutions), val(cluster_method)

    output:
    tuple val(sample), path("${sample}.doublets_removed.sct_clustered.rds"), emit: doublets_removed_rds
    tuple val(sample), path("${sample}.doublets_detected.rds"), emit: doublets_marked_rds
    tuple val(sample), path("qc_results"), emit: qc_results
    path "version.doubletfinder.txt", emit: version

    script:
    def mr = multiplet_rate == null ? '' : multiplet_rate
    assert default_res : 'Error: No default resolution supplied. Ensure each sample in in the samplesheet has a resolution under column name `res`.'
    def params_csv = 'param,value\n' +
        "multiplet_rate,${mr}\n" +
        "res,${default_res}\n" +
        "resolutions,${all_resolutions}\n" +
        "cluster_method,${cluster_method}\n"
    """
    # Create parameter samplesheet
    echo -e "${params_csv}" > params.csv

    doublet.R "${sample}" "${rds_path}" params.csv

    # Once doublets are removed, re-run SCTransform and clustering
    sct.R "${sample}" "${sample}.doublets_removed.rds" params.csv
    mv "${sample}.sct_clustered.rds" "${sample}.doublets_removed.sct_clustered.rds"

    # Get clustree version
    Rscript -e "cat(as.character(packageVersion('DoubletFinder')), '\\n')" > version.doubletfinder.txt
    """
}