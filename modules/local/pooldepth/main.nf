process POOLDEPTH {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/csvtk_mosdepth:64931f051298d5db':
        'community.wave.seqera.io/library/csvtk_mosdepth:e4d5a760c47164a4' }"

    input:
    tuple val(meta), path(bam), path(bai)
    val(pools) // list channel.of([1,2])
    val(window) // int - default should be 20
    
    output:
    tuple val(meta), path("${prefix}.depth.tsv"), emit: tsv
    tuple val("${task.process}"), val('pooldepth'), eval("csvtk version"), topic: versions, emit: versions_pooldepth
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir pool${pools[0]}
    mkdir pool${pools[1]}

    mosdepth pool${pools[0]}/${prefix} ${bam} --by ${window} --read-groups ${pools[0]}
    mosdepth pool${pools[1]}/${prefix} ${bam} --by ${window} --read-groups ${pools[1]}

    gunzip -c pool${pools[0]}/${prefix}.regions.bed.gz \\
    | csvtk mutate2 -H -t -n sample -e '"${prefix}"' \\
    | csvtk mutate2 -H -t -n pool -e '"${pools[0]}"' \\
    | csvtk add-header -H -t -n chrome,start,end,depth,sample,pool \\
    --out-file ${prefix}.${pools[0]}.depth.txt
    
    gunzip -c pool${pools[1]}/${prefix}.regions.bed.gz \\
    | csvtk mutate2 -H -t -n sample -e '"${prefix}"' \\
    | csvtk mutate2 -H -t -n pool -e '"${pools[1]}"' \\
    | csvtk add-header -H -t -n chrome,start,end,depth,sample,pool \\
    --out-file ${prefix}.${pools[1]}.depth.txt

    csvtk concat ${prefix}.${pools[0]}.depth.txt ${prefix}.${pools[1]}.depth.txt > ${prefix}.depth.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir pool${pools[0]}
    mkdir pool${pools[1]}

    touch ${prefix}.${pools[0]}.depth.txt
    touch ${prefix}.${pools[1]}.depth.txt
    touch ${prefix}.depth.tsv
    """
}
