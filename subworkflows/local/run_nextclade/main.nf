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

// mimics NEXTCLADE_DATASETGET version.yml
process grepTag {
    input:
    path dataset

    output:
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat <<-END_VERSIONS > versions.yml
    NEXTCLADE_DATASETGET:
        nextclade-tag: \$(grep "tag" $dataset/pathogen.json | sed -n 's/.*"tag": "\\([0-9-]\\+Z\\)".*/\\1/p')
    END_VERSIONS
    """
}

workflow RUN_NEXTCLADE {
    take: 
        ch_consensus
    
    main:
        nextclade_tag = params.nextclade_tag ?: ""
        ch_versions = channel.empty()

        if (file(params.nextclade).exists()) {
            nc_dataset = params.nextclade
            grepTag(nc_dataset)
            nc_tag = grepTag.out.versions
        } else if (params.nextclade instanceof String) {
            NEXTCLADE_DATASETGET(params.nextclade, nextclade_tag)
            nc_dataset = NEXTCLADE_DATASETGET.out.dataset
            nc_tag = NEXTCLADE_DATASETGET.out.versions
        }
        NEXTCLADE_RUN(ch_consensus, nc_dataset)
        ch_versions = ch_versions.mix(NEXTCLADE_RUN.out.versions)
        ch_versions = ch_versions.mix(nc_tag)

    emit:
        versions = ch_versions
        tsv = NEXTCLADE_RUN.out.tsv
}