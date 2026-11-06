process ARTIC_GUPPYPLEX {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/artic_gzip:9b6b316374421185'
        : 'artic/fieldbioinformatics:1.11.2'}"

    input:
    tuple val(meta), path(fastq_dir, stageAs: 'guppyplex_input/*')

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    tuple val("${task.process}"), val('artic'), eval('artic -v 2>&1 | sed "s/^.*artic //; s/ .*$//"'), emit: versions_artic, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Input may be a directory of FASTQs or FASTQ file(s) supplied directly. Inputs are staged into a
    // subdirectory so a file named '<prefix>.fastq(.gz)' can never collide with (and be overwritten by) the output.
    """
    if [ -d "${fastq_dir}" ]; then
        input_dir="${fastq_dir}"
    else
        input_dir="guppyplex_input"
        # guppyplex only collects '*.fastq*' files, so rename '.fq' / '.fq.gz' symlinks
        for f in guppyplex_input/*.fq guppyplex_input/*.fq.gz; do
            [ -e "\$f" ] || continue
            case "\$f" in
                *.fq.gz) mv "\$f" "\${f%.fq.gz}.fastq.gz" ;;
                *.fq)    mv "\$f" "\${f%.fq}.fastq" ;;
            esac
        done
    fi

    artic \\
        guppyplex \\
        ${args} \\
        --threads ${task.cpus} \\
        --directory \${input_dir} \\
        --output ${prefix}.fastq.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '' | gzip > ${prefix}.fastq.gz
    """
}
