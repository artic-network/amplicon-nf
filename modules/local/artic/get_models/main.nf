process ARTIC_GET_MODELS {
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/artic:1.6.4--pyhdfd78af_0'
        : 'biocontainers/artic:1.6.4--pyhdfd78af_0'}"

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
