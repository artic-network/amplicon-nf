/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FIND_CONCATENATE     } from '../../../modules/nf-core/find/concatenate'
include { NEXTCLADE_RUN        } from '../../../modules/nf-core/nextclade/run'
include { NEXTCLADE_DATASETGET } from '../../../modules/nf-core/nextclade/datasetget'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process grepTag {
    input:
    path dataset

    output:
    tuple val("${task.process}"), val('nextclade-tag'), eval("grep 'tag' ${dataset}/pathogen.json | sed -n 's/.*\"tag\": \"\\([0-9-]\\+Z\\)\".*/\\1/p'"), emit: versions_tag, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
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
    }
    else if (params.nextclade instanceof String) {
        NEXTCLADE_DATASETGET(params.nextclade, nextclade_tag)
        nc_dataset = NEXTCLADE_DATASETGET.out.dataset
    }
    NEXTCLADE_RUN(ch_consensus, nc_dataset)

    emit:
    versions = ch_versions
    tsv      = NEXTCLADE_RUN.out.tsv
}
