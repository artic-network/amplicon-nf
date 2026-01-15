/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CAT_CAT }                                 from '../../../modules/nf-core/cat/cat'
include { NEXTCLADE_RUN }                           from '../../../modules/nf-core/nextclade/run'
include { NEXTCLADE_DATASETGET }                    from '../../../modules/nf-core/nextclade/datasetget'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow RUN_NEXTCLADE {
    take: 
        ch_consensus
    
    main:
        nextclade_tag = params.nextclade_tag ?: ""
        NEXTCLADE_DATASETGET(params.nextclade, nextclade_tag)
        NEXTCLADE_RUN(ch_consensus, NEXTCLADE_DATASETGET.out.dataset)
    
    emit:
        versions = NEXTCLADE_RUN.out.versions
        tsv = NEXTCLADE_RUN.out.tsv
}