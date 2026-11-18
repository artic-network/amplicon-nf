process DEACON_FILTER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/deacon:0.19.1--h3edb6b3_0'
        : 'quay.io/biocontainers/deacon:0.19.1--h3edb6b3_0'}"

    input:
    // nanopore: fastq directory or single fastq file | illumina: [fastq_1, fastq_2]
    // Staged into a subdirectory so an input named like an output can never be overwritten
    tuple val(meta), path(reads, stageAs: 'input/*')
    path index

    output:
    tuple val(meta), path("*.dehosted.fastq.gz"), emit: fastq
    tuple val("${task.process}"), val('deacon'), eval('deacon --version | head -n1 | sed "s/deacon //"'), emit: versions_deacon, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Never add --rename: ONT read headers carry basecall_model_version_id, used to pick the Clair3 model
    if (meta.platform == 'nanopore') {
        """
        if [ -d ${reads} ]; then
            # Mirror guppyplex's input rules (top-level *.fastq*, skip *.temp) and stream everything through one deacon run
            # awful, will fix when deacon accepts gzipped input directly over stdin or directory input
            for f in ${reads}/*.fastq*; do
                case "\$f" in
                    *.temp) continue ;;
                    *.gz) gzip -dc "\$f" ;;
                    *) cat "\$f" ;;
                esac
            done | deacon \\
                filter \\
                --deplete \\
                --threads ${task.cpus} \\
                ${args} \\
                ${index} \\
                - \\
                --output ${prefix}.dehosted.fastq.gz
        else
            deacon \\
                filter \\
                --deplete \\
                --threads ${task.cpus} \\
                ${args} \\
                ${index} \\
                ${reads} \\
                --output ${prefix}.dehosted.fastq.gz
        fi
        """
    }
    else {
        """
        deacon \\
            filter \\
            --deplete \\
            --threads ${task.cpus} \\
            ${args} \\
            ${index} \\
            ${reads[0]} \\
            ${reads[1]} \\
            --output ${prefix}_R1.dehosted.fastq.gz \\
            --output2 ${prefix}_R2.dehosted.fastq.gz
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outputs = meta.platform == 'nanopore'
        ? "${prefix}.dehosted.fastq.gz"
        : "${prefix}_R1.dehosted.fastq.gz ${prefix}_R2.dehosted.fastq.gz"
    """
    for f in ${outputs}; do
        echo "" | gzip > \$f
    done
    """
}
