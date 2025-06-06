/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_amplicon-nf_pipeline'

include { ONT_ASSEMBLY           } from '../subworkflows/local/ont_assembly/main'
include { ILLUMINA_ASSEMBLY      } from '../subworkflows/local/illumina_assembly/main'

include { SAMTOOLS_DEPTH         } from '../modules/nf-core/samtools/depth/main'
include { SAMTOOLS_COVERAGE      } from '../modules/nf-core/samtools/coverage/main'

include { GENERATE_SAMPLE_REPORT } from '../modules/local/generate_sample_report/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow AMPLICON_NF {
    take:
    ch_samplesheet     // channel: samplesheet read in from --input
    ch_store_directory // channel: store directory read in from --store_dir

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_samplesheet
        .branch { meta, _fastq_directory, _fastq_1, _fastq_2 ->
            nanopore: meta.platform == "nanopore"
            illumina: meta.platform == "illumina"
        }
        .set { ch_input }

    //
    // Generate virus assemblies
    //
    ch_input.nanopore
        .map { meta, fastq_dir, _fastq_1, _fastq_2 ->
            [meta, fastq_dir]
        }
        .set { ch_nanopore_input }

    ONT_ASSEMBLY(
        ch_nanopore_input,
        ch_store_directory,
        ch_versions,
    )
    ch_versions = ch_versions.mix(ONT_ASSEMBLY.out.versions)

    ch_input.illumina
        .map { meta, _fastq_dir, fastq_1, fastq_2 ->
            [meta, fastq_1, fastq_2]
        }
        .set { ch_illumina_input }

    ILLUMINA_ASSEMBLY(
        ch_illumina_input,
        ch_store_directory,
        ch_versions,
    )
    ch_versions = ch_versions.mix(ILLUMINA_ASSEMBLY.out.versions)

    //
    // Generate report for each sample
    //
    ch_primertrimmed_bam = ONT_ASSEMBLY.out.primertrimmed_normalised_bam.mix(
        ILLUMINA_ASSEMBLY.out.primertrimmed_normalised_bam
    )

    ch_primer_scheme = ONT_ASSEMBLY.out.primer_scheme.mix(
        ILLUMINA_ASSEMBLY.out.primer_scheme
    )

    ch_amp_depth_tsv = ONT_ASSEMBLY.out.amplicon_depths.mix(
        ILLUMINA_ASSEMBLY.out.amplicon_depths
    )

    SAMTOOLS_COVERAGE(ch_primertrimmed_bam, [[:], []], [[:], []])
    ch_versions = ch_versions.mix(SAMTOOLS_COVERAGE.out.versions.first())

    ch_samtools_depth_input = ch_primertrimmed_bam.map { meta, bam, _bai ->
        [meta, bam]
    }

    SAMTOOLS_DEPTH(ch_samtools_depth_input, [[:], []])
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions.first())

    ch_sample_report_input = ch_primer_scheme
        .map { meta, bed, _ref -> [meta, bed] }
        .join(SAMTOOLS_DEPTH.out.tsv)
        .join(ch_amp_depth_tsv)
        .join(SAMTOOLS_COVERAGE.out.coverage)

    ch_report_template = file(
        "${projectDir}/assets/sample_report_template.html",
        checkIfExists: true
    )
    ch_artic_logo_svg = file(
        "${projectDir}/assets/artic-logo-small.svg",
        checkIfExists: true
    )
    ch_bootstrap_bundle_min_js = file(
        "${projectDir}/assets/bootstrap.bundle.min.js",
        checkIfExists: true
    )
    ch_bootstrap_bundle_min_css = file(
        "${projectDir}/assets/bootstrap.min.css",
        checkIfExists: true
    )
    ch_plotly_js = file(
        "${projectDir}/assets/plotly.min.js",
        checkIfExists: true
    )

    GENERATE_SAMPLE_REPORT(
        ch_sample_report_input,
        ch_report_template,
        ch_artic_logo_svg,
        ch_bootstrap_bundle_min_js,
        ch_bootstrap_bundle_min_css,
        ch_plotly_js,
    )
    ch_versions = ch_versions.mix(GENERATE_SAMPLE_REPORT.out.versions.first())

    ch_consensus_fasta = ONT_ASSEMBLY.out.consensus_fasta.mix(
        ILLUMINA_ASSEMBLY.out.consensus_fasta
    )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'amplicon-nf_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = Channel.fromPath(
        "${projectDir}/assets/multiqc_config.yml",
        checkIfExists: true
    )
    ch_multiqc_custom_config = params.multiqc_config
        ? Channel.fromPath(params.multiqc_config, checkIfExists: true)
        : Channel.empty()
    ch_multiqc_logo = params.multiqc_logo
        ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : Channel.empty()

    summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )
    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true,
        )
    )

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    emit:
    consensus_fasta = ch_consensus_fasta // channel: consensus FASTA files
    sample_report   = GENERATE_SAMPLE_REPORT.out.sample_report_html // channel: sample report files
    multiqc_report  = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions        = ch_versions // channel: software versions used in the workflow
}
