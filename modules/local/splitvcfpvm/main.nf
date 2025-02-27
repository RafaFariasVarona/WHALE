process SPLITVCFPVM {
    tag "$meta.id"
    label 'process_single'

    // Conda is not supported
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf)
    val n_vcf_variants_split

    output:
    tuple val(meta), path("*split*.vcf"), emit: vcfs
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
    grep "#" ${vcf} > ${prefix}-split00.vcf 

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
}
