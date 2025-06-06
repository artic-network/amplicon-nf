process ARTIC_ALIGN_TRIM {
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/artic:1.7.1--pyhdfd78af_0'
        : 'biocontainers/artic:1.7.1--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(sorted_bam), path(scheme_bed), path(scheme_ref)

    output:
    tuple val(meta), path("${prefix}.primertrimmed.sorted.bam"), emit: primertrimmed_bam
    tuple val(meta), path("${prefix}.amplicon_depths.tsv"), emit: amplicon_depths
    tuple val(meta), path("${prefix}.align_trim_report.tsv"), emit: alignreport
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"

    incorrect_pair_string = params.remove_incorrect_primer_pairs ? "--remove-incorrect-pairs" : ""
    endedness_string = meta.single_end ? "" : "--paired"

    """
    align_trim \\
        ${args} \\
        --normalise ${params.normalise_depth} \\
        --report ${prefix}.align_trim_report.tsv \\
        --amp-depth-report ${prefix}.amplicon_depths.tsv \\
        --trim-primers \\
        ${endedness_string} \\
        ${incorrect_pair_string} \\
        ${scheme_bed} \\
        < ${sorted_bam} \\
        > ${prefix}.primertrimmed.sam

    samtools sort -T ${prefix} ${prefix}.primertrimmed.sam -o ${prefix}.primertrimmed.sorted.bam


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(artic -v 2>&1 | sed 's/^.*artic //; s/ .*\$//')
    END_VERSIONS
    """
}
