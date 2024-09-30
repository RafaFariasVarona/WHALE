process SNIFFLES {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sniffles:2.3.3--pyhdfd78af_0' :
        'biocontainers/sniffles:2.3.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(gzi)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sniffles \\
        --input ${bam} \\
        --vcf ${prefix}_sniffles.vcf.gz \\
        --reference ${fasta} \\
        -t ${task.cpus} \\
        --sample-id ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sniffles: \$(sniffles --help 2>&1 | grep Version | sed 's/^.*Version //')
    END_VERSIONS
    """
}

