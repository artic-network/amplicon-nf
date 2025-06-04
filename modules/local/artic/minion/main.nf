process ARTIC_MINION {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/artic:1.6.4--pyhdfd78af_0'
        : 'biocontainers/artic:1.6.4--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(fastq), path(custom_scheme_directory)
    path store_directory

    output:
    tuple val(meta), path("${prefix}.sorted.bam"), path("${prefix}.sorted.bam.bai"), emit: sorted_bam
    tuple val(meta), path("${prefix}.trimmed.rg.sorted.bam"), path("${prefix}.trimmed.rg.sorted.bam.bai"), emit: normalised_bam
    tuple val(meta), path("${prefix}.primertrimmed.rg.sorted.bam"), path("${prefix}.primertrimmed.rg.sorted.bam.bai"), emit: primertrimmed_normalised_bam
    tuple val(meta), path("${prefix}.amplicon_depths.tsv"), emit: amplicon_depths
    tuple val(meta), path("${prefix}.consensus.fasta"), emit: fasta
    tuple val(meta), path("${prefix}.normalised.vcf.gz"), emit: vcf
    tuple val(meta), path("${prefix}.minion.log.txt"), emit: minion_log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"

    scheme_split = meta.scheme ? meta.scheme.split("/") : ["", "", ""]
    scheme_string = custom_scheme_directory ? "--bed ${custom_scheme_directory}/primer.bed --ref ${custom_scheme_directory}/reference.fasta" : "--scheme-name ${scheme_split[0]} --scheme-length ${scheme_split[1]} --scheme-version ${scheme_split[2]}"

    model_str = params.manual_clair3_model ? "--model ${params.manual_clair3_model}" : ""

    """
    artic \\
        minion \\
        ${args} \\
        --model-dir ${store_directory}/fieldbioinformatics-nf/clair3-models/ \\
        --scheme-directory ${store_directory}/fieldbioinformatics-nf/primer-schemes/ \\
        --threads ${task.cpus} \\
        --read-file ${fastq} \\
        ${model_str} \\
        ${scheme_string} \\
        ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(artic -v 2>&1 | sed 's/^.*artic //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.1.trimmed.rg.sorted.bam
    touch ${prefix}.1.trimmed.rg.sorted.bai
    touch ${prefix}.1.vcf
    touch ${prefix}.2.trimmed.rg.sorted.bam
    touch ${prefix}.2.trimmed.rg.sorted.bai
    touch ${prefix}.2.vcf

    touch ${prefix}.alignreport.csv
    touch ${prefix}.amplicon_depths.tsv

    touch ${prefix}.consensus.fasta
    touch ${prefix}.coverage_mask.txt
    touch ${prefix}.coverage_mask.txt.1.depths
    touch ${prefix}.coverage_mask.txt.2.depths

    touch ${prefix}.fail.vcf

    touch ${prefix}.merged.vcf
    echo "" | gzip > ${prefix}.merged.vcf.gz
    touch ${prefix}.merged.vcf.tbi

    touch ${prefix}.minion.log.txt

    echo "" | gzip > ${prefix}.normalised.vcf.gz
    touch ${prefix}.normalised.vcf.tbi

    touch ${prefix}.pass.vcf
    echo "" | gzip > ${prefix}.pass.vcf.gz
    touch ${prefix}.pass.vcf.gz.tbi

    touch ${prefix}.preconsensus.fasta
    touch ${prefix}.preconsensus.fasta.fai

    touch ${prefix}.primers.vcf
    touch ${prefix}.primersitereport.txt
    touch ${prefix}.primertrimmed.rg.sorted.bam
    touch ${prefix}.primertrimmed.rg.sorted.bam.bai

    touch ${prefix}.sorted.bam
    touch ${prefix}.sorted.bam.bai
    touch ${prefix}.trimmed.rg.sorted.bam
    touch ${prefix}.trimmed.rg.sorted.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(artic -v 2>&1 | sed 's/^.*artic //; s/ .*\$//')
    END_VERSIONS
    """
}
