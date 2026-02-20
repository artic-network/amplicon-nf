include { FASTCAT      } from '../../../modules/local/fastcat/main'
include { MOSDEPTH     } from '../../../modules/nf-core/mosdepth/main'

workflow AMPLICON_DEPTHS {

    take:
    ch_reads       // chnalle [meta, *fastq] 
    ch_bam_bai_bed // channel [meta, bam, bai, []]

    main:
    ch_versions = Channel.empty()
    
    FASTCAT(ch_reads)
    ch_versions = ch_versions.mix(FASTCAT.out.versions)
    println("finished fastcat")
    MOSDEPTH(ch_bam_bai_bed, [[],[]])
    println("finished mosdepth")


    emit:
    perfile  = FASTCAT.out.perfile    
    bed      = MOSDEPTH.out.per_base_bed
    versions = ch_versions
}
