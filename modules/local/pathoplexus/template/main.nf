process PATHOPLEXUS_TEMPLATE {
    tag "${meta.scheme ?: meta.custom_scheme_name}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/openpyxl_python:aa16c046cd10ee5c'
        : 'community.wave.seqera.io/library/openpyxl_python:7a42e0b9733f16d0'}"

    input:
    // jsons / fastas: PATHOPLEXUS_SAMPLE outputs for every sample of one primer scheme
    tuple val(meta), path(jsons, stageAs: 'samples/*'), path(fastas, stageAs: 'fastas/*')
    path template_config // assets/pathoplexus_template.toml

    output:
    tuple val(meta), path("${prefix}_pathoplexus_metadata.xlsx"), path("${prefix}_pathoplexus_metadata.tsv"), path("${prefix}_pathoplexus_sequences.fasta"), emit: upload
    // eval outputs are not allowed with python template scripts, so the script writes versions.yml
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: (meta.scheme ? meta.scheme.replace('/', '_') : meta.custom_scheme_name ? meta.custom_scheme_name.replace('/', '_') : 'unknown_scheme')
    template("pathoplexus_template.py")

    stub:
    prefix = task.ext.prefix ?: (meta.scheme ? meta.scheme.replace('/', '_') : meta.custom_scheme_name ? meta.custom_scheme_name.replace('/', '_') : 'unknown_scheme')
    """
    touch ${prefix}_pathoplexus_metadata.xlsx
    touch ${prefix}_pathoplexus_metadata.tsv
    touch ${prefix}_pathoplexus_sequences.fasta
    printf '"${task.process}":\\n  python: stub\\n  openpyxl: stub\\n' > versions.yml
    """
}
