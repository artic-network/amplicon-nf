process CSVTK_MUTATE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/csvtk:0.31.0--h9ee0642_0':
        'biocontainers/csvtk:0.31.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(csv, name: 'inputs/csv*/*')
    val in_format
    val out_format
    val subcommand

    output:
    tuple val(meta), path("${prefix}.${out_extension}"), emit: csv
    path "versions.yml"                                , emit: versions
    tuple val("${task.process}"), val('csvtk'), eval("csvtk version"), topic: versions, emit: versions_csvtk

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def delimiter = in_format == "tsv" ? "\t" : (in_format == "csv" ? "," : in_format)
    def out_delimiter = out_format == "tsv" ? "\t" : (out_format == "csv" ? "," : out_format)
    out_extension = out_format == "tsv" ? 'tsv' : 'csv'
    if (['mutate','mutate2','mutate3'].contains(subcommand)) {command = subcommand} else {error "Invalid function specified: ${subcommand}, must be one of: mutate, mutate2, mutate3"}

    """
    csvtk \\
        ${command} \\
        ${args} \\
        --num-cpus $task.cpus \\
        --delimiter "${delimiter}" \\
        --out-delimiter "${out_delimiter}" \\
        --out-file ${prefix}.${out_extension} \\
        ${csv}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvtk: \$(echo \$( csvtk version | sed -e "s/csvtk v//g" ))
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    out_extension = out_format == "tsv" ? 'tsv' : 'csv'
    """
    touch ${prefix}.${out_extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvtk: \$(echo \$( csvtk version | sed -e "s/csvtk v//g" ))
    END_VERSIONS
    """
}
