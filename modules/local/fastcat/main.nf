process FASTCAT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/fastcat:0.24.2--defceef4425784bd':
        'community.wave.seqera.io/library/fastcat:0.24.2--defceef4425784bd' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}.perfile.tsv"),             emit: perfile
    tuple val(meta), path("${prefix}.runids.tsv"),              emit: runids
    tuple val(meta), path("${prefix}.basecaller.tsv"),          emit: basecaller
    tuple val(meta), path("${prefix}.read.tsv"),                emit: read
    tuple val(meta), path("histograms/length.hist"),  emit: hist_length
    tuple val(meta), path("histograms/quality.hist"), emit: hist_quality
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
        --file ${prefix}.perfile.tsv \\
        --runids ${prefix}.runids.tsv \\
        --histograms histograms \\
        --basecallers ${prefix}.basecaller.tsv \\
        --read ${prefix}.read.tsv \\
        ${reads} > /dev/null 2>&1

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
