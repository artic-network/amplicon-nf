process GENERATE_SAMPLE_REPORT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/primalbedtools_jinja2_pandas_plotly:993c0aa5175f7683'
        : 'community.wave.seqera.io/library/jinja2_pandas_plotly_pip_primalbedtools:2b267b5fa768ec3d'}"

    input:
    tuple val(meta), path(bed), path(depth_tsv), path(amp_depth_tsv), path(coverage_report)
    path report_template
    path artic_logo_svg
    path bootstrap_bundle_min_js
    path bootstrap_bundle_min_css
    path plotly_js

    output:
    path "*_fieldbioinformatics-nf_report.html", emit: sample_report_html
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template("generate_sample_report.py")
}
