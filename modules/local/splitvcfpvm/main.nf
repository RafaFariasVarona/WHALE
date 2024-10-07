process SPLITVCFPVM {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf)
    val n_vcf_variants_split

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*split*.vcf"), emit: vcfs
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Hacemos el split
    grep -v "#" ${vcf} | split -a 3 -l ${n_vcf_variants_split} - split_vcf_


    n_files=\$(ls -l split_vcf_* | wc -l)
    let vcf_index=\${n_files}-1

    # Cambiamos el nombre a los archivos para que sean 01,02..
    counter=0
    for file in split_vcf_*; do
        mv "\$file" \$(printf "split_vcf_%02d" \$counter)
        counter=\$((\${counter} + 1))
    done

    # Guardamos la cabecera y creamos tantas cabeceras como archivos split tengamos.
    zgrep "#" ${vcf} > ${prefix}-split00.vcf 

    counter=1
    until [ \${counter} -gt \${vcf_index} ]
        do
        
        cp ${prefix}-split00.vcf ${prefix}-split0\${counter}.vcf
        counter=\$((\${counter} + 1))
    done

    # pegamos los split a la cabecera
    counter=0
    until [ \${counter} -gt \${vcf_index}  ]
        do
        echo \${counter}
        cat split_vcf_0\${counter} >> ${prefix}-split0\${counter}.vcf
        ##rm split_vcf_0\${counter}

        counter=\$((\${counter} + 1))
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        splitvcfpvm: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        splitvcfpvm: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
