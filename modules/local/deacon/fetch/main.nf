process DEACON_FETCH {
    tag "${index_name}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/deacon:0.20.0--h3edb6b3_0'
        : 'quay.io/biocontainers/deacon:0.20.0--h3edb6b3_0'}"

    input:
    path store_directory
    val index_name

    output:
    path store_directory, emit: store_directory
    tuple val("${task.process}"), val('deacon'), eval('deacon --version | head -n1 | sed "s/deacon //"'), emit: versions_deacon, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def index_dir = "${store_directory}/amplicon-nf/deacon/"
    // Short-term: deacon (openssl-probe 0.2) doesn't look in conda's ssl/ dir for CA certs
    def ca_cert_file = workflow.containerEngine ? '/usr/local/ssl/cacert.pem' : '${CONDA_PREFIX}/ssl/cacert.pem'
    """
    export SSL_CERT_FILE=${ca_cert_file}

    mkdir -p ${index_dir}

    # deacon (>= 0.20) skips the download if the index already exists (add --force via ext.args to
    # replace it) and never leaves a truncated download in place of a complete index
    deacon \\
        index fetch \\
        ${args} \\
        ${index_name} \\
        --output ${index_dir}/${index_name}.idx
    """

    stub:
    def index_dir = "${store_directory}/amplicon-nf/deacon/"
    """
    mkdir -p ${index_dir}
    touch ${index_dir}/${index_name}.idx
    """
}
