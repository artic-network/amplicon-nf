process FASTCAT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/fastcat:0.24.2--defceef4425784bd':
        'community.wave.seqera.io/library/fastcat:0.24.2--defceef4425784bd' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${prefix}.perfile.tsv"),             optional:true, emit: summary
    tuple val(meta), path("${prefix}.runid.tsv"),               optional:true, emit: runid
    tuple val(meta), path("${prefix}.basecaller.tsv"),          optional:true, emit: basecaller
    tuple val(meta), path("${prefix}.read.tsv"),                optional:true, emit: read
    tuple val(meta), path("${prefix}/histograms/length.hist"),  optional:true, emit: hist_length
    tuple val(meta), path("${prefix}/histograms/quality.hist"), optional:true, emit: hist_quality
    tuple val("${task.process}"), val('fastcat'), eval("fastcat --version"), topic: versions, emit: versions_fastcat
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    fastcat \\
        $args \\
        --threads $task.cpus \\
        --sample ${prefix} \\
        $fastq 2&>1 /dev/null

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastcat: \$(fastcat --version 2>&1))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args
    
    touch ${prefix}.perfile.tsv
    touch ${prefix}.runid.tsv
    touch ${prefix}.basecaller.tsv
    touch ${prefix}.read.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastcat: \$(fastcat --version 2>&1)
    END_VERSIONS
    """
}
