process SVIM {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/svim:2.0.0--pyhdfd78af_0':
        'biocontainers/svim:2.0.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(gzi)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    svim alignment \\
        . \\
        ${bam} \\
        ${fasta} \\
        --sample ${prefix} \\
        ${args}

    mv variants.vcf ${prefix}_svim.vcf
    gzip ${prefix}_svim.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svim: \$(svim --version | sed 's/svim //')
    END_VERSIONS
    """
}
