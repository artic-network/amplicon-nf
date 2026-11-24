process ARTIC_GET_SCHEME {
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0ff516afa828b5ba2c95c5c88c870938d146ae2961f6c6b26ab8d00b21567c94/data'
        : 'artic/fieldbioinformatics:1.12.0'}"

    input:
    tuple val(meta), path(fastq_1), path(fastq_2)
    path store_directory

    output:
    tuple val(meta), path(fastq_1), path(fastq_2), path("primer.bed"), path("reference.fasta"), emit: reads_and_scheme
    tuple val("${task.process}"), val('artic'), eval('artic -v 2>&1 | sed "s/^.*artic //; s/ .*$//"'), emit: versions_artic, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    scheme_split = meta.scheme.split("/")

    """
    artic_get_scheme \\
        --scheme-directory ${store_directory}/amplicon-nf/primer-schemes/ \\
        --scheme-name ${scheme_split[0]} \\
        --scheme-length ${scheme_split[1]} \\
        --scheme-version ${scheme_split[2]} \\
        --read-file ${fastq_1}
    """
}
