include { BWAMEM2_MEM                                        } from '../../../modules/nf-core/bwamem2/mem/main'
include { BWAMEM2_INDEX                                      } from '../../../modules/nf-core/bwamem2/index/main'
include { SAMTOOLS_SORT                                      } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX                                     } from '../../../modules/nf-core/samtools/index/main'
include { TABIX_BGZIPTABIX                                   } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { TRIMMOMATIC                                        } from '../../../modules/nf-core/trimmomatic/main'
include { SAMTOOLS_FAIDX                                     } from '../../../modules/nf-core/samtools/faidx/main'
include { FREEBAYES                                          } from '../../../modules/nf-core/freebayes/main'
include { BCFTOOLS_NORM                                      } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_CONSENSUS as BCFTOOLS_CONSENSUS_AMBIGUOUS } from '../../../modules/nf-core/bcftools/consensus/main'
include { BCFTOOLS_CONSENSUS as BCFTOOLS_CONSENSUS_FIXED     } from '../../../modules/nf-core/bcftools/consensus/main'
include { BCFTOOLS_VIEW                                      } from '../../../modules/nf-core/bcftools/view/main'

include { ARTIC_GET_SCHEME                                   } from '../../../modules/local/artic/get_scheme/main'
include { ARTIC_ALIGNTRIM                                    } from '../../../modules/nf-core/artic/aligntrim/main'
include { PROCESS_GVCF                                       } from '../../../modules/local/process_gvcf/main'
include { SEQKIT_REPLACE_IUPAC                               } from '../../../modules/local/seqkit/replace_iupac/main'


workflow ILLUMINA_ASSEMBLY {
    take:
    ch_input
    ch_store_directory

    main:

    ch_branched_input = ch_input.branch { meta, _fastq_1, _fastq_2 ->
        remote_scheme: meta.scheme
        custom_scheme: meta.custom_scheme
    }

    ch_custom_scheme_input = ch_branched_input.custom_scheme.map { meta, fastq_1, fastq_2 ->
        [meta, fastq_1, fastq_2, file("${meta.custom_scheme.toUriString()}/primer.bed", checkIfExists: true), file("${meta.custom_scheme.toUriString()}/reference.fasta", checkIfExists: true)]
    }

    ARTIC_GET_SCHEME(
        ch_branched_input.remote_scheme,
        ch_store_directory,
    )

    ch_reads_and_scheme = ARTIC_GET_SCHEME.out.reads_and_scheme.mix(ch_custom_scheme_input)

    ch_trimmomatic_input = ch_reads_and_scheme.map { meta, fastq_1, fastq_2, _scheme_bed, _scheme_ref ->
        [meta, [fastq_1, fastq_2]]
    }

    TRIMMOMATIC(ch_trimmomatic_input)

    ch_trimmed_fastq = TRIMMOMATIC.out.trimmed_reads
        .map { meta, trimmed_fastq ->
            [meta, trimmed_fastq[0], trimmed_fastq[1]]
        }
        .join(
            ch_reads_and_scheme.map { meta, _fastq_1, _fastq_2, scheme_bed, scheme_ref ->
                [meta, scheme_bed, scheme_ref]
            }
        )

    // No need to index a ref for each sample, just once per scheme
    ch_refs_only = ch_reads_and_scheme
        .map { meta, _fastq_1, _fastq_2, _scheme_bed, scheme_ref ->
            [meta.subMap("scheme", "custom_scheme", "custom_scheme_name"), scheme_ref]
        }
        .groupTuple()
        .map { meta, refs -> [meta, refs[0]] }

    // Replace IUPAC ambiguity codes in reference FASTA — freebayes does not support them
    SEQKIT_REPLACE_IUPAC(ch_refs_only)

    SEQKIT_REPLACE_IUPAC.out.iupac_replaced.subscribe { _flag ->
        log.warn("IUPAC ambiguity codes were detected in one or more reference sequences and have been replaced with N. This is because Freebayes does not support IUPAC ambiguity codes. Please check the output FASTA files to ensure that this is acceptable for your analysis.")
    }

    ch_refs_only = SEQKIT_REPLACE_IUPAC.out.fasta

    // Propagate cleaned reference into ch_trimmed_fastq so all downstream tools (bwa-mem2, freebayes, bcftools) receive it
    ch_trimmed_fastq = ch_trimmed_fastq
        .map { meta, fastq_1, fastq_2, scheme_bed, _old_ref ->
            [meta.subMap("scheme", "custom_scheme", "custom_scheme_name"), meta, fastq_1, fastq_2, scheme_bed]
        }
        .combine(ch_refs_only, by: 0)
        .map { _scheme_meta, meta, fastq_1, fastq_2, scheme_bed, cleaned_ref ->
            [meta, fastq_1, fastq_2, scheme_bed, cleaned_ref]
        }

    BWAMEM2_INDEX(
        ch_refs_only
    )

    // Join the index with the per-sample metadata
    ch_rejoined_bwamem2_indices = BWAMEM2_INDEX.out.index
        .combine(
            ch_reads_and_scheme.map { meta, _fastq_1, _fastq_2, _scheme_bed, _scheme_ref ->
                [meta.subMap("scheme", "custom_scheme", "custom_scheme_name"), meta]
            },
            by: 0
        )
        .map { _scheme_meta, bwamem2_index, meta ->
            [meta, bwamem2_index]
        }

    ch_bwamem2_mem_input = ch_trimmed_fastq
        .join(ch_rejoined_bwamem2_indices)
        .multiMap { meta, fastq_1, fastq_2, _scheme_bed, scheme_ref, scheme_ref_index ->
            reads: [meta, [fastq_1, fastq_2]]
            ref_index: [meta, scheme_ref_index]
            ref_fasta: [meta, scheme_ref]
        }

    BWAMEM2_MEM(
        ch_bwamem2_mem_input.reads,
        ch_bwamem2_mem_input.ref_index,
        ch_bwamem2_mem_input.ref_fasta,
        true,
    )

    ch_sorted_bam = BWAMEM2_MEM.out.bam.join(
        ch_trimmed_fastq.map { meta, _fastq_1, _fastq_2, scheme_bed, scheme_ref ->
            [meta, scheme_bed, scheme_ref]
        }
    )

    ch_aligntrim_input = ch_sorted_bam.map { meta, sorted_bam, scheme_bed, _scheme_ref ->
        [meta, sorted_bam, scheme_bed, params.normalise_depth ?: []]
    }

    // Sort the output BAMfile
    ARTIC_ALIGNTRIM(ch_aligntrim_input, true)

    SAMTOOLS_INDEX(ARTIC_ALIGNTRIM.out.primertrimmed_bam)

    ch_samtools_faidx_input = ch_refs_only.map { meta, scheme_ref ->
        [meta, scheme_ref, []]
    }

    SAMTOOLS_FAIDX(ch_samtools_faidx_input, false)

    // Join the fai with the per-sample metadata
    ch_rejoined_ref_fais = SAMTOOLS_FAIDX.out.fai
        .combine(
            ch_reads_and_scheme.map { meta, _fastq_1, _fastq_2, _scheme_bed, _scheme_ref ->
                [meta.subMap("scheme", "custom_scheme", "custom_scheme_name"), meta]
            },
            by: 0
        )
        .map { _scheme_meta, scheme_ref_fai, meta ->
            [meta, scheme_ref_fai]
        }

    ch_freebayes_input = ARTIC_ALIGNTRIM.out.primertrimmed_bam
        .join(SAMTOOLS_INDEX.out.index)
        .join(
            ch_sorted_bam.map { meta, _sorted_bam, scheme_bed, scheme_ref ->
                [meta, scheme_bed, scheme_ref]
            }
        )
        .join(ch_rejoined_ref_fais)
        .multiMap { meta, primertrimmed_bam, primertrimmed_bam_bai, _scheme_bed, scheme_ref, scheme_ref_fai ->
            sam_input: [meta, primertrimmed_bam, primertrimmed_bam_bai, [], [], []]
            ref_input: [meta, scheme_ref]
            ref_fai_input: [meta, scheme_ref_fai]
        }

    // LOTS of stuff to do in modules.conf for this to work
    FREEBAYES(ch_freebayes_input.sam_input, ch_freebayes_input.ref_input, ch_freebayes_input.ref_fai_input, [[:], []], [[:], []], [[:], []])

    PROCESS_GVCF(
        FREEBAYES.out.vcf
    )

    ch_bcftools_norm_input = PROCESS_GVCF.out.consensus_vcf
        .join(ch_freebayes_input.ref_input)
        .multiMap { meta, vcf, scheme_ref ->
            vcf_input: [meta, vcf, []]
            ref_input: [meta, scheme_ref]
        }

    BCFTOOLS_NORM(ch_bcftools_norm_input.vcf_input, ch_bcftools_norm_input.ref_input)

    ch_vartype = Channel.of("fixed", "ambiguous")

    ch_bcftools_view_input = BCFTOOLS_NORM.out.vcf
        .join(BCFTOOLS_NORM.out.tbi)
        .combine(ch_vartype)
        .map { meta, vcf, vcf_index, vartype_str ->
            [meta + [vartype: vartype_str], vcf, vcf_index]
        }


    BCFTOOLS_VIEW(ch_bcftools_view_input, [], [], [])

    ch_bcftools_consensus_reference = ch_bcftools_norm_input.ref_input
        .combine(ch_vartype)
        .map { meta, scheme_ref, vartype_str ->
            [meta + [vartype: vartype_str], scheme_ref]
        }

    ch_bcftools_consensus_depth_mask = PROCESS_GVCF.out.depth_mask
        .combine(ch_vartype)
        .map { meta, depth_mask, vartype_str ->
            [meta + [vartype: vartype_str], depth_mask]
        }

    ch_bcftools_consensus_input = BCFTOOLS_VIEW.out.vcf
        .join(BCFTOOLS_VIEW.out.tbi)
        .join(ch_bcftools_consensus_reference)
        .join(ch_bcftools_consensus_depth_mask)
        .branch { meta, _vcf, _vcf_index, _reference, _depth_mask ->
            fixed: meta.vartype == "fixed"
            ambiguous: meta.vartype == "ambiguous"
        }

    // Don't use the depth mask for ambiguous variants preconsensus generation
    ch_bcftools_consensus_input_ambiguous = ch_bcftools_consensus_input.ambiguous.map { meta, vcf, vcf_index, reference, _depth_mask ->
        [meta, vcf, vcf_index, reference, []]
    }

    BCFTOOLS_CONSENSUS_AMBIGUOUS(ch_bcftools_consensus_input_ambiguous)

    ch_preconsensus_fasta = BCFTOOLS_CONSENSUS_AMBIGUOUS.out.fasta.map { meta, fasta ->
        [meta - meta.subMap("vartype"), fasta]
    }

    ch_bcftools_consensus_input_fixed = ch_bcftools_consensus_input.fixed
        .map { meta, vcf, vcf_index, reference, depth_mask ->
            [meta - meta.subMap("vartype"), vcf, vcf_index, reference, depth_mask]
        }
        .join(ch_preconsensus_fasta)
        .map { meta, vcf, vcf_index, _scheme_ref, depth_mask, preconsensus ->
            [meta, vcf, vcf_index, preconsensus, depth_mask]
        }

    BCFTOOLS_CONSENSUS_FIXED(ch_bcftools_consensus_input_fixed)

    // Join the primertrimmed bam with its index
    ch_primertrimmed_bam = ARTIC_ALIGNTRIM.out.primertrimmed_bam.join(
        SAMTOOLS_INDEX.out.index
    )

    ch_primer_scheme = ch_reads_and_scheme.map { meta, _fastq_1, _fastq_2, scheme_bed, scheme_ref ->
        [meta, scheme_bed, scheme_ref]
    }

    emit:
    consensus_fasta              = BCFTOOLS_CONSENSUS_FIXED.out.fasta
    amplicon_depths              = ARTIC_ALIGNTRIM.out.amp_depth_report
    sorted_bam                   = BWAMEM2_MEM.out.bam
    primertrimmed_normalised_bam = ch_primertrimmed_bam
    primer_scheme                = ch_primer_scheme
}
