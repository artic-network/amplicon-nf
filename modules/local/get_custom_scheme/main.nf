process GET_CUSTOM_SCHEME {
    tag "${meta.id}"
    label 'process_single'

    //  This container is used elsewhere in the pipeline, so we use it here too since it just needs any container
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pysam:0.23.0--py312h47d5410_0'
        : 'biocontainers/pysam:0.23.0--py312h47d5410_0'}"

    input:
    tuple val(meta), path(fastq_1), path(fastq_2), path(custom_scheme_directory)

    output:
    tuple val(meta), path(fastq_1), path(fastq_2), path("${custom_scheme_directory}/primer.bed"), path("${custom_scheme_directory}/reference.fasta"), emit: reads_and_scheme

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    """
}
