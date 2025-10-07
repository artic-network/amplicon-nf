process GENERATE_RUN_REPORT {
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/pip_jinja2_biopython_numpy_pruned:62f9c8b9b1f39c59'
        : 'community.wave.seqera.io/library/pip_jinja2_biopython_numpy_pruned:fbb9affc06c37839'}"

    input:
    tuple val(meta), path(bed), path(depth_tsvs, stageAs: "depth_tsvs/*"), path(amp_depth_tsvs, stageAs: "amplicon_depth_tsvs/*"), path(coverage_tsvs, stageAs: "coverage_tsvs/*"), path(msas, stageAs: "msas/*"), path(samplesheet_csv)
    path report_template
    path artic_logo_svg
    path bootstrap_bundle_min_js
    path bootstrap_bundle_min_css
    path plotly_js

    output:
    path "*_amplicon-nf_run-report.html", emit: run_report_html
    path "*_qc_results.tsv", emit: qc_results_tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template("generate_run_report.py")
}
