process ARTIC_GUPPYPLEX {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/artic:1.10.3--aa87beb2acab0fa4'
        : 'artic/fieldbioinformatics:1.10.3'}"

    input:
    tuple val(meta), path(fastq_dir)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    tuple val("${task.process}"), val('artic'), eval('artic -v 2>&1 | sed "s/^.*artic //; s/ .*$//"'), emit: versions_artic, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    artic \\
        guppyplex \\
        ${args} \\
        --threads ${task.cpus} \\
        --directory ${fastq_dir} \\
        --output ${prefix}.fastq

    pigz -p ${task.cpus} *.fastq
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '' | gzip > ${prefix}.fastq.gz
    """
}
