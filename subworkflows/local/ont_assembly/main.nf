include { ARTIC_GUPPYPLEX  } from '../../../modules/nf-core/artic/guppyplex/main'
include { ARTIC_MINION     } from '../../../modules/local/artic/minion/main'
include { ARTIC_GET_MODELS } from '../../../modules/local/artic/get_models/main'

workflow ONT_ASSEMBLY {
    take:
    ch_input
    ch_store_directory

    main:

    ch_input
        .count()
        .filter { it > 0 }
        .combine(Channel.value(ch_store_directory))
        .map { _count, store_dir -> store_dir }
        .first()
        .set { ch_gated_store_directory }

    ARTIC_GET_MODELS(ch_gated_store_directory)

    ARTIC_GUPPYPLEX(
        ch_input
    )

    ch_guppyplexed_fastq = ARTIC_GUPPYPLEX.out.fastq

    ch_minion_input = ch_guppyplexed_fastq.map { meta, fastq_dir ->
        [meta, fastq_dir, meta.custom_scheme ? file(meta.custom_scheme) : []]
    }

    ARTIC_MINION(
        ch_minion_input,
        ARTIC_GET_MODELS.out.store_directory,
    )

    emit:
    consensus_fasta              = ARTIC_MINION.out.fasta
    amplicon_depths              = ARTIC_MINION.out.amplicon_depths
    pre_normalisation_coverage   = ARTIC_MINION.out.pre_normalisation_coverage_report
    post_normalisation_coverage  = ARTIC_MINION.out.post_normalisation_coverage_report
    sorted_bam                   = ARTIC_MINION.out.sorted_bam
    primertrimmed_normalised_bam = ARTIC_MINION.out.primertrimmed_normalised_bam
    primer_scheme                = ARTIC_MINION.out.primer_scheme
}
