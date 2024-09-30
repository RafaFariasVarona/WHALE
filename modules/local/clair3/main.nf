process CLAIR3 {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/clair3:1.0.8--py39hf5e1c6e_1':
        'biocontainers/clair3:1.0.8--py39hf5e1c6e_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(gzi)

    output:
    tuple val(meta), path("*clair3.vcf.gz"), emit: vcf
    tuple val(meta), path("*clair3.vcf.gz.tbi"), emit: vcf_tbi
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    run_clair3.sh \\
        --bam_fn=${bam} \\
        --ref_fn=${fasta} \\
        --threads=${task.cpus} \\
        ${args} \\
        --output=. \\
        --sample_name=${prefix}.C3
    
    mv merge_output.vcf.gz ${prefix}_clair3.vcf.gz
    mv merge_output.vcf.gz.tbi ${prefix}_clair3.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(clair3 --version |& sed '1!d ; s/clair3 //')
    END_VERSIONS
    """
}
