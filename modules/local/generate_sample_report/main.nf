process GENERATE_SAMPLE_REPORT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/primalbedtools_biopython_jinja2_numpy_pruned:2ad9362062fb41eb'
        : 'community.wave.seqera.io/library/primalbedtools_biopython_jinja2_numpy_pruned:d3e3819de3f6e323'}"

    input:
    tuple val(meta), path(bed), path(depth_tsv), path(amp_depth_tsv), path(coverage_report)
    path report_template
    path artic_logo_svg
    path bootstrap_bundle_min_js
    path bootstrap_bundle_min_css
    path plotly_js

    output:
    path "*_amplicon-nf_sample-report.html", emit: sample_report_html
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template("generate_sample_report.py")
}
