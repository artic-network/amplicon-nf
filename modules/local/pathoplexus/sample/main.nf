process PATHOPLEXUS_SAMPLE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/openpyxl_python:aa16c046cd10ee5c'
        : 'community.wave.seqera.io/library/openpyxl_python:7a42e0b9733f16d0'}"

    input:
    // consensus: reheadered consensus FASTA | bed: per-sample primer.bed | qc_tsv: per-sample QC report TSV
    // reads: reads after the dehosting step (dehosted, or raw if dehosting was skipped)
    tuple val(meta), path(consensus, stageAs: 'consensus/*'), path(bed, stageAs: 'scheme/*'), path(qc_tsv, stageAs: 'qc/*'), path(reads, stageAs: 'reads/*')
    val pipeline_version  // e.g. 0.12.0-g44f4001, see getPipelineVersion()
    val nondefault_params // free text for the `comment` field, see nonDefaultParamsText()

    output:
    tuple val(meta), path("${meta.id}.pathoplexus.json"), path("${meta.id}.pathoplexus.fasta"), emit: sample
    // Header-stripped copies of the dehosted reads of included samples (nothing for excluded samples)
    tuple val(meta), path("raw_reads/*.dehosted.fastq.gz"), emit: raw_reads, optional: true
    // eval outputs are not allowed with python template scripts, so the script writes versions.yml
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // JSON-encoded so the text is always a valid Python string literal in the template, whatever characters it holds
    nondefault_params_json = groovy.json.JsonOutput.toJson(nondefault_params)
    template("pathoplexus_sample.py")

    stub:
    """
    echo '{"id": "${meta.id}", "excluded": false, "fasta_ids": ["${meta.id}"], "row": {"id": "${meta.id}"}, "raw_read_files": []}' > ${meta.id}.pathoplexus.json
    printf '>${meta.id}\\nN\\n' > ${meta.id}.pathoplexus.fasta
    printf '"${task.process}":\\n  python: stub\\n' > versions.yml
    """
}
