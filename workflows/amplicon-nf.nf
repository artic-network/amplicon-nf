/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                                   } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                    } from '../subworkflows/local/utils_nfcore_amplicon-nf_pipeline'

include { ONT_ASSEMBLY                              } from '../subworkflows/local/ont_assembly/main'
include { ILLUMINA_ASSEMBLY                         } from '../subworkflows/local/illumina_assembly/main'

include { SAMTOOLS_DEPTH                            } from '../modules/nf-core/samtools/depth/main'
include { SAMTOOLS_COVERAGE                         } from '../modules/nf-core/samtools/coverage/main'
include { SEQKIT_REPLACE as SEQKIT_REPLACE_ONT      } from '../modules/nf-core/seqkit/replace/main'
include { SEQKIT_REPLACE as SEQKIT_REPLACE_ILLUMINA } from '../modules/nf-core/seqkit/replace/main'
include { MAFFT_ALIGN                               } from '../modules/nf-core/mafft/align/main'
include { SEQKIT_GREP as SEQKIT_GREP_FASTAS         } from '../modules/nf-core/seqkit/grep/main'
include { SEQKIT_GREP as SEQKIT_GREP_REFS           } from '../modules/nf-core/seqkit/grep/main'
include { CAT_CAT                                   } from '../modules/nf-core/cat/cat/main'

include { GENERATE_SAMPLE_REPORT                    } from '../modules/local/generate_sample_report/main'
include { GENERATE_RUN_REPORT                       } from '../modules/local/generate_run_report/main'

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

    ch_input = ch_samplesheet.branch { meta, _fastq_directory, _fastq_1, _fastq_2 ->
        nanopore: meta.platform == "nanopore"
        illumina: meta.platform == "illumina"
    }

    if (params.read_directory) {
        // 
        // Illumina file fuzzy matching from read_directory
        //
        read_directory = file(params.read_directory, checkIfExists: true)

        ch_illumina_missing_files = ch_input.illumina
            .filter { _meta, _fastq_dir, fastq_1, fastq_2 ->
                !fastq_1 || !fastq_2
            }
            .map { meta, _fastq_dir, _fastq_1, _fastq_2 ->
                [meta.id, meta]
            }

        // Find Illumina file names from the fastq_directory, then get rid of the sample suffix (e.g. _S1) from the file names so it'll match the sample names in the samplesheet.
        ch_file_pairs = Channel.fromFilePairs(
                "${read_directory}/*_R{1,2}*.fastq.gz",
                size: 2,
                type: "file",
                followLinks: true,
                maxDepth: 1,
            )
            .map { common_key, file_pair ->
                [common_key.replaceAll("_L00.\$", "").replaceAll("_S.\$", ""), file_pair[0], file_pair[1]]
            }

        ch_fuzzy_matched_illumina = ch_illumina_missing_files
            .join(ch_file_pairs, failOnDuplicate: true)
            .map { _sample_name, meta, fastq_1, fastq_2 ->
                [meta, fastq_1, fastq_2]
            }

        // 
        // ONT file fuzzy matching from read_directory
        //
        ch_nanopore_fuzzy_match = ch_input.nanopore
            .filter { meta, fastq_dir, _fastq_1, _fastq_2 ->
                !fastq_dir && meta.barcode
            }
            .map { meta, _fastq_dir, _fastq_1, _fastq_2 ->
                [meta, file("${read_directory}/${meta.barcode}/", checkIfExists: true)]
            }

        // Mix in the fuzzy matched Nanopore files with the well behaved samplesheet entries.
        ch_nanopore_input = ch_input.nanopore
            .filter { meta, fastq_dir, _fastq_1, _fastq_2 ->
                fastq_dir && !meta.barcode
            }
            .map { meta, fastq_dir, _fastq_1, _fastq_2 ->
                [meta, fastq_dir]
            }
            .mix(ch_nanopore_fuzzy_match)

        // Mix in the fuzzy matched Illumina files with the well behaved samplesheet entries.
        ch_illumina_input = ch_input.illumina
            .filter { _meta, _fastq_dir, fastq_1, fastq_2 ->
                fastq_1 && fastq_2
            }
            .map { meta, _fastq_dir, fastq_1, fastq_2 ->
                [meta, fastq_1, fastq_2]
            }
            .mix(ch_fuzzy_matched_illumina)
    }
    else {
        ch_nanopore_input = ch_input.nanopore
            .filter { meta, fastq_dir, _fastq_1, _fastq_2 ->
                fastq_dir && !meta.barcode
            }
            .map { meta, fastq_dir, _fastq_1, _fastq_2 ->
                [meta, fastq_dir]
            }

        ch_illumina_input = ch_input.illumina
            .filter { _meta, _fastq_dir, fastq_1, fastq_2 ->
                fastq_1 && fastq_2
            }
            .map { meta, _fastq_dir, fastq_1, fastq_2 ->
                [meta, fastq_1, fastq_2]
            }
    }

    // If there's no work to do, produce an informative error message stating the possible causes.
    ch_nanopore_input
        .mix(ch_illumina_input)
        .ifEmpty {
            error(
                """
                No valid input found. Please check the following:
                
                1) If you are using implicit (fuzzy) input matching, ensure you are providing the `--read_directory` parameter and that the directory structure / file naming conventions are compatible with the expected patterns as described in https://github.com/artic-network/amplicon-nf/blob/main/docs/usage.md, this is crucial for successful file matching. This is the most common cause of this error being raised.

                2) The samplesheet provided via --input is not empty.
                """
            )
        }
        .filter { it != null }


    //
    // Generate virus assemblies
    //    

    ONT_ASSEMBLY(
        ch_nanopore_input,
        ch_store_directory,
        ch_versions,
    )
    ch_versions = ch_versions.mix(ONT_ASSEMBLY.out.versions)

    SEQKIT_REPLACE_ONT(ONT_ASSEMBLY.out.consensus_fasta)
    ch_versions = ch_versions.mix(SEQKIT_REPLACE_ONT.out.versions.first())

    ILLUMINA_ASSEMBLY(
        ch_illumina_input,
        ch_store_directory,
        ch_versions,
    )
    ch_versions = ch_versions.mix(ILLUMINA_ASSEMBLY.out.versions)

    SEQKIT_REPLACE_ILLUMINA(ILLUMINA_ASSEMBLY.out.consensus_fasta)
    ch_versions = ch_versions.mix(SEQKIT_REPLACE_ILLUMINA.out.versions.first())

    ch_reheadered_consensus_fasta = SEQKIT_REPLACE_ILLUMINA.out.fastx.mix(
        SEQKIT_REPLACE_ONT.out.fastx
    )

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

    sample_report_template = file(
        "${projectDir}/assets/sample_report_template.html",
        checkIfExists: true
    )
    run_report_template = file(
        "${projectDir}/assets/run_report_template.html",
        checkIfExists: true
    )
    artic_logo_svg = file(
        "${projectDir}/assets/artic-logo-small.svg",
        checkIfExists: true
    )
    bootstrap_bundle_min_js = file(
        "${projectDir}/assets/bootstrap.bundle.min.js",
        checkIfExists: true
    )
    bootstrap_bundle_min_css = file(
        "${projectDir}/assets/bootstrap.min.css",
        checkIfExists: true
    )
    plotly_js = file(
        "${projectDir}/assets/plotly.min.js",
        checkIfExists: true
    )

    GENERATE_SAMPLE_REPORT(
        ch_sample_report_input,
        sample_report_template,
        artic_logo_svg,
        bootstrap_bundle_min_js,
        bootstrap_bundle_min_css,
        plotly_js,
    )
    ch_versions = ch_versions.mix(GENERATE_SAMPLE_REPORT.out.versions.first())

    ch_bed_by_scheme = ch_primer_scheme
        .map { meta, bed, _ref ->
            [
                meta.subMap("scheme", "custom_scheme", "custom_scheme_name"),
                bed,
            ]
        }
        .unique()

    if (params.primer_mismatch_plot) {
        ch_chroms = ch_bed_by_scheme
            .splitCsv(elem: 1, header: false, sep: "\t", strip: true)
            .filter { _meta, bed_row ->
                !bed_row[0].toString().startsWith("#")
            }
            .map { meta, bed_row -> [meta, bed_row[0]] }
            .unique()

        ch_consensus_by_chrom = ch_chroms
            .combine(
                ch_reheadered_consensus_fasta.map { meta, fasta -> [meta.subMap("scheme", "custom_scheme", "custom_scheme_name"), fasta] },
                by: 0
            )
            .map { meta, chrom, fasta ->
                [
                    meta + [chrom: chrom] + [id: chrom],
                    fasta,
                ]
            }
            .groupTuple()

        CAT_CAT(ch_consensus_by_chrom)
        ch_versions = ch_versions.mix(CAT_CAT.out.versions.first())

        SEQKIT_GREP_FASTAS(CAT_CAT.out.file_out, [])
        ch_versions = ch_versions.mix(SEQKIT_GREP_FASTAS.out.versions.first())

        ch_refs_per_chrom = ch_chroms
            .combine(ch_primer_scheme.map { meta, _bed, ref -> [meta.subMap("scheme", "custom_scheme", "custom_scheme_name"), ref] }, by: 0)
            .map { meta, chrom, ref -> [meta + [chrom: chrom, id: chrom], ref] }

        SEQKIT_GREP_REFS(ch_refs_per_chrom, [])
        ch_versions = ch_versions.mix(SEQKIT_GREP_REFS.out.versions.first())

        ch_mafft_align_input = SEQKIT_GREP_FASTAS.out.filter
            .join(SEQKIT_GREP_REFS.out.filter)
            .multiMap { meta, fastas, reference ->
                fastas: [meta, fastas]
                reference: [meta, reference]
            }

        MAFFT_ALIGN(ch_mafft_align_input.reference, [[:], []], ch_mafft_align_input.fastas, [[:], []], [[:], []], [[:], []], false)
        ch_versions = ch_versions.mix(MAFFT_ALIGN.out.versions.first())

        ch_msas_by_scheme = MAFFT_ALIGN.out.fas
            .map { meta, msa ->
                [
                    meta.subMap("scheme", "custom_scheme", "custom_scheme_name"),
                    msa,
                ]
            }
            .groupTuple()
    }

    ch_amp_depth_tsvs_by_scheme = ch_amp_depth_tsv
        .map { meta, tsv -> [meta.subMap("scheme", "custom_scheme", "custom_scheme_name"), tsv] }
        .groupTuple()

    ch_depth_tsvs_by_scheme = SAMTOOLS_DEPTH.out.tsv
        .map { meta, tsv -> [meta.subMap("scheme", "custom_scheme", "custom_scheme_name"), tsv] }
        .groupTuple()

    ch_coverage_tsvs_by_scheme = SAMTOOLS_COVERAGE.out.coverage
        .map { meta, tsv -> [meta.subMap("scheme", "custom_scheme", "custom_scheme_name"), tsv] }
        .groupTuple()

    samplesheet_csv = file("${params.input}", checkIfExists: true)

    // optional reporting outputs
    ch_msas_opt       = params.primer_mismatch_plot ? ch_msas_by_scheme : channel.value([])
    ch_nextclade_opt  = params.nextclade ? ch_nextclade_tsv : channel.value([])

    ch_run_report_input = ch_bed_by_scheme
        // defaults
        .join(ch_depth_tsvs_by_scheme)
        .join(ch_amp_depth_tsvs_by_scheme)
        .join(ch_coverage_tsvs_by_scheme)
        // optional
        .join(ch_msas_opt)
        .join(ch_nextclade_opt)
        .map { meta, bed, depth, amp, cov, msas, nextclade ->
            [
                meta,
                bed,
                depth,
                amp,
                cov,
                msas ?: [],
                nextclade ?: [],
                samplesheet_csv
            ]
        }

    // if (params.primer_mismatch_plot) {
    //     ch_run_report_input = ch_bed_by_scheme
    //         .join(ch_depth_tsvs_by_scheme)
    //         .join(ch_amp_depth_tsvs_by_scheme)
    //         .join(ch_coverage_tsvs_by_scheme)
    //         .join(ch_msas_by_scheme)
    //         .map { meta, bed, depth_tsvs, amp_depth_tsvs, coverage_tsvs, msas ->
    //             [
    //                 meta,
    //                 bed,
    //                 depth_tsvs,
    //                 amp_depth_tsvs,
    //                 coverage_tsvs,
    //                 msas,
    //                 samplesheet_csv,
    //             ]
    //         }
    // }
    // else {
    //     ch_run_report_input = ch_bed_by_scheme
    //         .join(ch_depth_tsvs_by_scheme)
    //         .join(ch_amp_depth_tsvs_by_scheme)
    //         .join(ch_coverage_tsvs_by_scheme)
    //         .map { meta, bed, depth_tsvs, amp_depth_tsvs, coverage_tsvs ->
    //             [meta, bed, depth_tsvs, amp_depth_tsvs, coverage_tsvs, [], samplesheet_csv]
    //         }
    // }

    GENERATE_RUN_REPORT(
        ch_run_report_input,
        run_report_template,
        artic_logo_svg,
        bootstrap_bundle_min_js,
        bootstrap_bundle_min_css,
        plotly_js,
    )
    ch_versions = ch_versions.mix(GENERATE_RUN_REPORT.out.versions.first())

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
    multiqc_report  = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions        = ch_versions // channel: software versions used in the workflow    
    consensus_fasta = ch_reheadered_consensus_fasta // channel: consensus FASTA files
    sample_report   = GENERATE_SAMPLE_REPORT.out.sample_report_html // channel: sample report files
}
