process DEEPVARIANT {
    tag "$meta.id"
    label 'process_high'

    //Conda is not supported at the moment
    container "docker://quay.io/nf-core/deepvariant:1.6.1"

    input:
    tuple val(meta), path(input), path(index)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(gzi)

    output:
    tuple val(meta), path("${prefix}_deepvariant.vcf.gz")      ,  emit: vcf
    tuple val(meta), path("${prefix}_deepvariant.vcf.gz.tbi")  ,  emit: vcf_tbi
    tuple val(meta), path("${prefix}_deepvariant.g.vcf.gz")    ,  emit: gvcf
    tuple val(meta), path("${prefix}_deepvariant.g.vcf.gz.tbi"),  emit: gvcf_tbi
    path "versions.yml"                            ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DEEPVARIANT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    /opt/deepvariant/bin/run_deepvariant \\
        --ref=${fasta} \\
        --reads=${input} \\
        --output_vcf=${prefix}_deepvariant.vcf.gz \\
        --output_gvcf=${prefix}_deepvariant.g.vcf.gz \\
        --intermediate_results_dir=tmp \\
        --num_shards=${task.cpus} \\
        ${args} \\
        --sample_name=${prefix}.DV

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant: \$(echo \$(/opt/deepvariant/bin/run_deepvariant --version) | sed 's/^.*version //; s/ .*\$//' )
    END_VERSIONS
    """
}
