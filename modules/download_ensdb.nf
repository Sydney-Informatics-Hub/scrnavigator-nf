process DOWNLOAD_ENSDB {
    publishDir "${params.outdir}/ensdb", mode: 'copy'
    container "sydneyinformaticshub/scrnavigator-nf-annotate"

    input:
    val(species)
    val(ref_version)

    output:
    path("EnsDb_*.sqlite"), emit: ensdb

    script:
    """
    mkdir -p .cache
    export XDG_CACHE_HOME="\${PWD}/.cache"
    
    download_ensdb.R "${species}" "${ref_version}" "\${XDG_CACHE_HOME}"
    """
}
