include { FASTCAT } from '../../../modules/local/fastcat/main'
include { CSVTK_CONCAT as CONCAT_SUMMARY } from '../../../modules/nf-core/csvtk/concat/main'
include { CSVTK_CONCAT as CONCAT_DEPTH } from '../../../modules/nf-core/csvtk/concat/main'
include { POOLDEPTH } from "../../../modules/local/pooldepth/main"

workflow AMPLICON_DEPTHS {
    take:
    ch_reads // chnalle [meta, *fastq] 
    ch_bam_bai // channel [meta, bam, bai]
    pools // channel.of([1,2])
    window // val(20)

    main:
    ch_versions = channel.empty()

    FASTCAT(ch_reads)
    perfile_collected = FASTCAT.out.perfile
        .collect { it[1] }
        .map { [[id: "all"], it] }
    CONCAT_SUMMARY(perfile_collected, "tsv", "tsv")

    POOLDEPTH(ch_bam_bai, pools, window)
    depth_collected = POOLDEPTH.out.tsv
        .collect { it[1] }
        .map { [[id: "all"], it] }
    CONCAT_DEPTH(depth_collected, "tsv", "tsv")

    ch_versions = ch_versions
        .mix(FASTCAT.out.versions)
        .mix(CONCAT_SUMMARY.out.versions)
        .mix(POOLDEPTH.out.versions)

    emit:
    summary = CONCAT_SUMMARY.out.csv
    bed = CONCAT_DEPTH.out.csv
    versions = ch_versions
}