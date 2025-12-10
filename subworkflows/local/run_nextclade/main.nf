/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CAT_CAT }                                 from '../../../modules/nf-core/cat/cat'
include { NEXTCLADE_DATASETGET }                    from '../../../modules/nf-core/nextclade/datasetget'
include { NEXTCLADE_RUN }                           from '../../../modules/nf-core/nextclade/run'
include { SEQKIT_REPLACE as SEQKIT_REPLACE_NC }     from "../../../modules/nf-core/seqkit/replace"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow RUN_NEXTCLADE {
    take: 
        ch_consensus // channel from amplicon-nf.nf 
    
    main:
        nextclade_tag = params.nextclade_tag ?: ""
        
        NEXTCLADE_DATASETGET (
            params.nextclade,
            nextclade_tag
        )

        ch_all_consensus_fasta = ch_consensus
            .map { _meta, fasta -> fasta }
            .collectFile(name: 'all_consensus.fasta')
            .map { multi_fasta ->[[id: 'all_consensus.fasta'], multi_fasta]}
        
        SEQKIT_REPLACE_NC(ch_all_consensus_fasta)

        NEXTCLADE_RUN (
            SEQKIT_REPLACE_NC.out.fastx,
            NEXTCLADE_DATASETGET.out.dataset
        )
    
    emit:
        versions = NEXTCLADE_RUN.out.versions
        nextclade_tsv = NEXTCLADE_RUN.out.tsv.map {_meta, tsv -> tsv}
}