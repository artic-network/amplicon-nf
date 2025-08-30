process ARTIC_MINION {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/artic:1.8.1--ec21cd6688fda876'
        : 'community.wave.seqera.io/library/artic:1.8.1--8ca66e901fbb9970'}"

    input:
    tuple val(meta), path(fastq), path(custom_scheme_directory)
    path store_directory

    output:
    tuple val(meta), path("${prefix}.sorted.bam"), path("${prefix}.sorted.bam.bai"), emit: sorted_bam
    tuple val(meta), path("${prefix}.primertrimmed.rg.sorted.bam"), path("${prefix}.primertrimmed.rg.sorted.bam.bai"), emit: primertrimmed_normalised_bam
    tuple val(meta), path("${prefix}.amplicon_depths.tsv"), emit: amplicon_depths
    tuple val(meta), path("${prefix}.consensus.fasta"), emit: fasta
    tuple val(meta), path("${prefix}.normalised.vcf.gz"), emit: vcf
    tuple val(meta), path("primer.bed"), path("reference.fasta"), emit: primer_scheme
    tuple val(meta), path("${prefix}.minion.log.txt"), emit: minion_log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"

    scheme_split = meta.scheme ? meta.scheme.split("/") : ["", "", ""]
    scheme_string = custom_scheme_directory ? "--bed ${custom_scheme_directory}/primer.bed --ref ${custom_scheme_directory}/reference.fasta" : "--scheme-name ${scheme_split[0]} --scheme-length ${scheme_split[1]} --scheme-version ${scheme_split[2]}"
    scheme_copy_string = custom_scheme_directory ? "cp ${custom_scheme_directory}/primer.bed primer.bed && cp ${custom_scheme_directory}/reference.fasta reference.fasta" : "cp ${store_directory}/amplicon-nf/primer-schemes/${scheme_split[0]}/${scheme_split[1]}/${scheme_split[2]}/primer.bed primer.bed && cp ${store_directory}/amplicon-nf/primer-schemes/${scheme_split[0]}/${scheme_split[1]}/${scheme_split[2]}/reference.fasta reference.fasta"

    """
    artic \\
        minion \\
        ${args} \\
        --model-dir ${store_directory}/amplicon-nf/clair3-models/ \\
        --scheme-directory ${store_directory}/amplicon-nf/primer-schemes/ \\
        --threads ${task.cpus} \\
        --read-file ${fastq} \\
        ${scheme_string} \\
        ${prefix}

    ${scheme_copy_string}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(artic -v 2>&1 | sed 's/^.*artic //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.1.primertrimmed.rg.sorted.bam
    touch ${prefix}.1.primertrimmed.rg.sorted.bam.bai
    touch ${prefix}.1.vcf
    touch ${prefix}.2.primertrimmed.rg.sorted.bam
    touch ${prefix}.2.primertrimmed.rg.sorted.bam.bai
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
