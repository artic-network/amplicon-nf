process ARTIC_GET_MODELS {
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/artic:1.10.1--194b4fac7363edea'
        : 'artic/fieldbioinformatics:1.10.1'}"

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
