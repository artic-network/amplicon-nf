include { ARTIC_GUPPYPLEX  } from '../../../modules/nf-core/artic/guppyplex/main'
include { ARTIC_MINION     } from '../../../modules/local/artic/minion/main'
include { ARTIC_GET_MODELS } from '../../../modules/local/artic/get_models/main'

workflow ONT_ASSEMBLY {
    take:
    ch_input
    ch_store_directory
    ch_versions

    main:

    ARTIC_GET_MODELS(ch_store_directory)
    ch_versions = ch_versions.mix(ARTIC_GET_MODELS.out.versions)

    ARTIC_GUPPYPLEX(
        ch_input
    )
    ch_versions = ch_versions.mix(ARTIC_GUPPYPLEX.out.versions.first())
    ARTIC_GUPPYPLEX.out.fastq.set { ch_guppyplexed_fastq }

    ch_guppyplexed_fastq
        .map { meta, fastq_dir ->
            [meta, fastq_dir, meta.custom_scheme ? file(meta.custom_scheme) : []]
        }
        .set { ch_minion_input }

    ARTIC_MINION(
        ch_minion_input,
        ARTIC_GET_MODELS.out.store_directory,
    )
    ch_versions = ch_versions.mix(ARTIC_MINION.out.versions.first())

    emit:
    consensus_fasta              = ARTIC_MINION.out.fasta
    amplicon_depths              = ARTIC_MINION.out.amplicon_depths
    sorted_bam                   = ARTIC_MINION.out.sorted_bam
    primertrimmed_normalised_bam = ARTIC_MINION.out.primertrimmed_normalised_bam
    primer_scheme                = ARTIC_MINION.out.primer_scheme
    versions                     = ch_versions
}
