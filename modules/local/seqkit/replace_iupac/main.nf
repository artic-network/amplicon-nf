process SEQKIT_REPLACE_IUPAC {
    tag "${meta.scheme ?: meta.custom_scheme_name ?: meta.toString()}"
    label 'process_low'

    conda "bioconda::seqkit=2.9.0"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0'
        : 'quay.io/biocontainers/seqkit:2.9.0--h9ee0642_0'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fasta"),             emit: fasta
    path "iupac_replaced.flag", optional: true,   emit: iupac_replaced
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/^.*v//'"), emit: versions_seqkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if grep -v '^>' "${fasta}" | grep -qiE '[RYWSKMBDHV]'; then
        seqkit \\
            replace \\
            -s \\
            -p '[RYWSKMBDHVrywskmbdhv]' \\
            -r 'N' \\
            --threads ${task.cpus} \\
            -i ${fasta} \\
            -o ${prefix}.fasta
        touch iupac_replaced.flag
    else
        cp ${fasta} ${prefix}.fasta
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" > ${prefix}.fasta
    """
}
