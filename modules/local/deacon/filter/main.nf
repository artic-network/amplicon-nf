process DEACON_FILTER {
    tag "${meta.id}"
    label 'process_low'

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
            # Stream every top-level FASTQ file through one deacon run. Only names ending .fastq / .fq
            # (optionally .gz, any case) are read, so sidecar files such as *.fastq.gz.deacon.json, *.temp,
            # hidden files and sub-directories are skipped, as are empty files.
            shopt -s nullglob nocasematch
            gz=()
            plain=()
            for f in ${reads}/*; do
                name="\${f##*/}"
                [ -f "\$f" ] && [ -s "\$f" ] || continue
                case "\$name" in
                    .*) continue ;;
                    *.fastq|*.fq|*.fastq.gz|*.fq.gz) ;;
                    *) continue ;;
                esac
                # Classify by the gzip magic bytes (1f 8b) rather than the name, in case the two disagree
                if [ "\$(head -c 2 "\$f" | od -An -tx1 | tr -d ' \\n')" = "1f8b" ]; then
                    gz+=("\$f")
                else
                    plain+=("\$f")
                fi
            done
            shopt -u nocasematch

            if [ \$(( \${#gz[@]} + \${#plain[@]} )) -eq 0 ]; then
                echo "ERROR: no non-empty FASTQ files (.fastq, .fq, .fastq.gz, .fq.gz) found in '${reads}'" >&2
                exit 1
            fi
            echo "Dehosting \${#gz[@]} gzipped and \${#plain[@]} uncompressed FASTQ file(s) from '${reads}'" >&2

            # deacon (>= 0.20) detects gzip on stdin and reads concatenated .gz files as one stream, so the
            # files are passed through compressed. One stream can't switch between gzip and plain text, so a
            # directory mixing the two is decompressed instead. Empty arrays are never expanded (bash -u).
            if [ \${#plain[@]} -eq 0 ]; then
                cat "\${gz[@]}"
            elif [ \${#gz[@]} -eq 0 ]; then
                cat "\${plain[@]}"
            else
                gzip -dc "\${gz[@]}"
                cat "\${plain[@]}"
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
