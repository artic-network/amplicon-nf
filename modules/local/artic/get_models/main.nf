process ARTIC_GET_MODELS {
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0ff516afa828b5ba2c95c5c88c870938d146ae2961f6c6b26ab8d00b21567c94/data'
        : 'artic/fieldbioinformatics:1.12.0'}"

    input:
    path store_directory

    output:
    path store_directory, emit: store_directory
    tuple val("${task.process}"), val('artic'), eval('artic -v 2>&1 | sed "s/^.*artic //; s/ .*$//"'), emit: versions_artic, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    artic_get_models \\
        --model-dir ${store_directory}/amplicon-nf/clair3-models/
    """
}
