process GET_VCF_FIELDS {
    tag "$meta.id"
    label 'process_single'

    // Conda is not supported
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${prefix}_GT.txt"), emit: gt
    tuple val(meta), path("${prefix}_DP_mean.txt"), emit: dp_mean
    tuple val(meta), path("${prefix}_AD_mean.txt"), emit: ad_mean
    tuple val(meta), path("${prefix}_VAF.txt"), emit: vaf
    tuple val(meta), path("${prefix}_caller_order.txt"), emit: caller_order
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    bcftools query -f '[\\t%GT]\\n' ${vcf} | sed 's/\\t//1' | sed 's/\\.\\/\\.//g' | sed 's/|/\\//g' | sed 's/1\\/0/0\\/1/g' > ${prefix}_GT.txt
	bcftools query -f '[\\t%DP]\\n' ${vcf} | sed 's/\\t//1' | sed 's/\\.//g' > ${prefix}_DP.txt
	bcftools query -f '[\\t%AD{1}]\\n' ${vcf} | sed 's/\\t//1' | sed 's/\\.//g' > ${prefix}_VD.txt
	bcftools query -f '[\\t%AD{0}]\\n' ${vcf} | sed 's/\\t//1' | sed 's/\\.//g' > ${prefix}_RD.txt

	awk -v OFMT=%.0f '{sum = 0; for (i = 1; i <= NF; i++) sum += \$i; sum /= NF; print sum}' ${prefix}_DP.txt > ${prefix}_DP_mean.txt
    awk -v OFMT=%.0f '{sum = 0; count = 0; for (i = 1; i <= NF; i++) if (\$i != 0) {sum += \$i; count++} if (count > 0) sum /= count; print sum}' ${prefix}_VD.txt > ${prefix}_VD_mean.txt
	awk -v OFMT=%.0f '{sum = 0; count = 0; for (i = 1; i <= NF; i++) if (\$i != 0) {sum += \$i; count++} if (count > 0) sum /= count; print sum}' ${prefix}_RD.txt > ${prefix}_RD_mean.txt
	
    paste -d, ${prefix}_RD_mean.txt ${prefix}_VD_mean.txt > ${prefix}_AD_mean.txt

    paste ${prefix}_VD_mean.txt ${prefix}_DP_mean.txt | awk -v OFMT=%.2f '{print(\$1/\$2)}' > ${prefix}_VAF.txt

    gzip -dc ${vcf} | grep '^#CHROM' | awk -F'\\t' '{for (i=1; i<=NF; i++) if (\$i ~ /\\./) print substr(\$i, index(\$i, ".")+1)}' > ${prefix}_caller_order.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
// Get genotype (GT), read depth (DP), variant allele depth (VD) and reference allele depth (RD)
// Calculate the average DP (DP_mean), VD (VD_mean) and RD (RD_mean)
// For VD_mean and RD_mean exclude 0 values from the average calculation
// Merge RD_mean.txt and VD_mean.txt with a comma delimiter and write the result to AD_mean.txt
// Calculate variant allele frequency (VAF)
// Get the variant callers in order of appearance 
