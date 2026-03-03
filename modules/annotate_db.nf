process ANNOTATE_DATABASE {
    publishDir "${params.outdir}/annotation/db", mode: 'copy'
    // ext input_size: { new InputFileSizes(rds_path) + new InputFileSizes(annotation_db) }
    ext input_size: { rds_path.size() + ( annotation_db ? annotation_db.size() : 0 ) }
    memory { 8.GB + task.ext.input_size.B * 10 }

    input:
    tuple val(cohort_name), path(rds_path), val(species), val(min_cells_for_annotation), path(annotation_db)

    output:
    tuple val(cohort_name), path("${cohort_name}.annotated.database.rds"), emit: annotated_rds
    tuple val(cohort_name), path("available_annotations.txt"), emit: available_annotations
    tuple val(cohort_name), path("qc_results"), emit: qc_results

    script:
    def annotation_csv = 'param,value\n' +
        "species,${species}\n" +
        "min_cells_for_annotation,${min_cells_for_annotation}\n" +
        "annotation_db,${annotation_db}\n"
    """
    # Create annotation samplesheet
    echo -e "${annotation_csv}" > annotation.csv

    # Create a cache dir
    mkdir -p .cache
    export XDG_CACHE_HOME="\${PWD}/.cache"

    annotate_db.R "${cohort_name}" "${rds_path}" annotation.csv
    """
}