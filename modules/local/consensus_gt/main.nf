process CONSENSUS_GT {
    tag "$meta.id"
    label 'process_single'

    // Conda is not supported
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1':
        'biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(GT)
    tuple val(meta), path(caller_order)

    output:
    tuple val(meta), path("${prefix}_GT_consensus.txt"), emit: gt_consensus
    tuple val(meta), path("${prefix}_GT_discordances.txt"), emit: gt_discordances
    
    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env python3
    
    import pandas as pd

    df = pd.read_table("${GT}", sep="\\t", header=None)
    
    order = pd.read_table("${caller_order}", sep="\\t", header=None)
    order["pos"] = range(0, order.shape[0])

    df_new = pd.concat([df[order[order[0] == "C3"]["pos"]]] * 3 + [df[order[order[0] == "DV"]["pos"]]] * 4 + [df[order[order[0] == "NC"]["pos"]]] * 2 , axis=1, ignore_index=True)
    
    df_new.mode(axis=1).to_csv("${prefix}_GT_consensus.txt", sep='\\t', header=False, index=False)
    
    discordance_list=list(map(lambda x: len(set(filter(lambda y: y == y , x))) - 1,df.values))
    with open("${prefix}_GT_discordances.txt", 'w') as fp:
        for item in discordance_list:
            fp.write("%s\\n" % item)
    """
}
// Python script (from PARROT-FJD) to find the consensus genotype between callers