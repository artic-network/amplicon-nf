process DEACON_FILTER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/deacon:0.20.0--h3edb6b3_0'
        : 'quay.io/biocontainers/deacon:0.20.0--h3edb6b3_0'}"

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
            files=()
            n_gz=0
            for f in ${reads}/*.fastq*; do
                case "\$f" in
                    *.temp) ;;
                    *.gz) files+=("\$f"); n_gz=\$((n_gz + 1)) ;;
                    *) files+=("\$f") ;;
                esac
            done
            # deacon (>= 0.20) detects gzip on stdin and reads concatenated .gz files as one stream, so the
            # files are passed through compressed. One stream can't switch between gzip and plain text, so a
            # directory mixing the two is decompressed file by file instead.
            if [ \$n_gz -eq 0 ] || [ \$n_gz -eq \${#files[@]} ]; then
                cat "\${files[@]}"
            else
                for f in "\${files[@]}"; do
                    case "\$f" in
                        *.gz) gzip -dc "\$f" ;;
                        *) cat "\$f" ;;
                    esac
                done
            fi | deacon \\
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
