process SAMPLE_INFO {
    tag "$meta.id"
    label 'process_single'

    // Conda is not supported
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(gt_consensus), path(ad_mean), path(dp_mean), path(vaf), path(sf), path(gt_discordances)

    output:
    tuple val(meta), path("${prefix}_sample_info.txt"), emit: sample_info
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    bcftools query -f '[:%GT:%DP:%AD{1}]\\n' ${vcf} > ${prefix}_format_SF.txt
	
    paste -d ":" ${gt_consensus} ${ad_mean} ${dp_mean} ${vaf} ${sf} ${gt_discordances} > ${prefix}_format_join.txt
	
    paste ${prefix}_format_join.txt ${prefix}_format_SF.txt | tr -d '\\t' > ${prefix}_sample_info.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
// Create SAMPLE column information