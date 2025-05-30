process PROCESS_GVCF {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pysam:0.23.0--py312h47d5410_0'
        : 'biocontainers/pysam:0.23.0--py312h47d5410_0'}"

    input:
    tuple val(meta), path(gvcf)

    output:
    tuple val(meta), path("variants.vcf"), emit: variants_vcf
    tuple val(meta), path("consensus.vcf"), emit: consensus_vcf
    tuple val(meta), path("mask.txt"), emit: depth_mask

    when:
    task.ext.when == null || task.ext.when

    script:
    template("process_gvcf.py")
}
