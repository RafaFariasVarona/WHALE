process HEADER_VARIANTS_VCF {
    tag "$meta.id"
    label 'process_single'

    // Conda is not supported
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(caller_order), path(format_vcf), path(sample_info)

    output:
    tuple val(meta), path("${prefix}_final.vcf"), emit: final_vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    bcftools view -h ${vcf} | grep "##" > ${prefix}_header.vcf

    echo "##FORMAT=<ID=SF,Number=1,Type=String,Description=\\"Software\\">" >> ${prefix}_header.vcf
	echo "##FORMAT=<ID=GD,Number=1,Type=String,Description=\\"Genotype discordances. 0 same genotype and 1 different genotype\\">" >> ${prefix}_header.vcf

	echo "##FORMAT=<ID="\$(awk 'NR==1 {print}' ${caller_order})"_GT,Number=1,Type=String,Description=\\""\$(awk 'NR==1 {print}' ${caller_order})" genotype\\">" >> ${prefix}_header.vcf
	echo "##FORMAT=<ID="\$(awk 'NR==1 {print}' ${caller_order})"_DP,Number=1,Type=String,Description=\\""\$(awk 'NR==1 {print}' ${caller_order})" read depth\\">" >> ${prefix}_header.vcf
	echo "##FORMAT=<ID="\$(awk 'NR==1 {print}' ${caller_order})"_VD,Number=1,Type=String,Description=\\""\$(awk 'NR==1 {print}' ${caller_order})" variant allele depth\\">" >> ${prefix}_header.vcf

	echo "##FORMAT=<ID="\$(awk 'NR==2 {print}' ${caller_order})"_GT,Number=1,Type=String,Description=\\""\$(awk 'NR==2 {print}' ${caller_order})" genotype\\">" >> ${prefix}_header.vcf
	echo "##FORMAT=<ID="\$(awk 'NR==2 {print}' ${caller_order})"_DP,Number=1,Type=String,Description=\\""\$(awk 'NR==2 {print}' ${caller_order})" read depth\\">" >> ${prefix}_header.vcf
	echo "##FORMAT=<ID="\$(awk 'NR==2 {print}' ${caller_order})"_VD,Number=1,Type=String,Description=\\""\$(awk 'NR==2 {print}' ${caller_order})" variant allele depth\\">" >> ${prefix}_header.vcf

	echo "##FORMAT=<ID="\$(awk 'NR==3 {print}' ${caller_order})"_GT,Number=1,Type=String,Description=\\""\$(awk 'NR==3 {print}' ${caller_order})" genotype\\">" >> ${prefix}_header.vcf
	echo "##FORMAT=<ID="\$(awk 'NR==3 {print}' ${caller_order})"_DP,Number=1,Type=String,Description=\\""\$(awk 'NR==3 {print}' ${caller_order})" read depth\\">" >> ${prefix}_header.vcf
	echo "##FORMAT=<ID="\$(awk 'NR==3 {print}' ${caller_order})"_VD,Number=1,Type=String,Description=\\""\$(awk 'NR==3 {print}' ${caller_order})" variant allele depth\\">" >> ${prefix}_header.vcf

    sed 's/\\(DV\\)/DeepVariant/2' ${prefix}_header.vcf > ${prefix}_header1.vcf
    sed 's/\\(C3\\)/Clair3/2' ${prefix}_header1.vcf >> ${prefix}_header2.vcf
    sed 's/\\(NC\\)/NanoCaller/2' ${prefix}_header2.vcf >> ${prefix}_final.vcf

    bcftools view -h ${vcf} | tail -n 1 | cut -f 1-10 | sed "s/."\$(awk 'NR==1 {print}' ${caller_order})"//" >> ${prefix}_final.vcf

    paste -d "\\t" ${format_vcf} ${sample_info} >> ${prefix}_final.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
// Add new information to the header of the final vcf and paste the new vcf info (from CHROM to FORMAT) and SAMPLE info