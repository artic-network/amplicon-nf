process ARTIC_GET_MODELS {
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/artic:1.9.0--4bfb8149af4d3e92'
        : 'community.wave.seqera.io/library/artic:1.9.0--ed3ab66c9589cea3'}"

    input:
    path store_directory

    output:
    path store_directory, emit: store_directory
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    artic_get_models \\
        --model-dir ${store_directory}/amplicon-nf/clair3-models/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(artic -v 2>&1 | sed 's/^.*artic //; s/ .*\$//')
    END_VERSIONS
    """
}
