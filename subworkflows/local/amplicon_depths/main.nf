include { FASTCAT } from '../../../modules/local/fastcat/main'
include { MOSDEPTH } from '../../../modules/nf-core/mosdepth/main'
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main'
include { CAT_CAT } from "../../../modules/nf-core/cat/cat/main"

workflow AMPLICON_DEPTHS {
    take:
    ch_reads // chnalle [meta, *fastq] 
    ch_bam_bai_bed // channel [meta, bam, bai, []]

    main:
    ch_versions = channel.empty()

    FASTCAT(ch_reads)
    perfile_collected = FASTCAT.out.perfile
        .collect { it[1] }
        .map { [[id: "all"], it] }
    CSVTK_CONCAT(perfile_collected, "tsv", "tsv")

    MOSDEPTH(ch_bam_bai_bed, [[], []])
    beds_collected = MOSDEPTH.out.per_base_bed
        .collect { it[1] }
        .map { [[id: "all"], it] }
    CAT_CAT(beds_collected)

    ch_versions = ch_versions
        .mix(FASTCAT.out.versions)
        .mix(CAT_CAT.out.versions)

    emit:
    summary = CSVTK_CONCAT.out.csv
    bed = CAT_CAT.out.file_out
    versions = ch_versions
}
