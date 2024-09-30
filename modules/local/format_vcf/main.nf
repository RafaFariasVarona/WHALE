process FORMAT_VCF {
    tag "$meta.id"
    label 'process_single'

    // Conda is not supported
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(caller_order)

    output:
    tuple val(meta), path("${prefix}_format_vcf.txt"), emit: format_vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    bcftools view -H ${vcf} | cut -f 1-8 | sed 's/\$/\\tGT:AD:DP:VAF:SF:GD:'"\$(awk 'NR==1 {print}' ${caller_order})"'_GT:'"\$(awk 'NR==1 {print}' ${caller_order})"'_DP:'"\$(awk 'NR==1 {print}' ${caller_order})"'_VD:'"\$(awk 'NR==2 {print}' ${caller_order})"'_GT:'"\$(awk 'NR==2 {print}' ${caller_order})"'_DP:'"\$(awk 'NR==2 {print}' ${caller_order})"'_VD:'"\$(awk 'NR==3 {print}' ${caller_order})"'_GT:'"\$(awk 'NR==3 {print}' ${caller_order})"'_DP:'"\$(awk 'NR==3 {print}' ${caller_order})"'_VD/' > ${prefix}_format_vcf.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
// Take CHROM, POS, ID, REF, ALT, QUAL, FILTER and INFO	columns and add a new FORMAT column = "GT:AD:DP:VAF:SF:GD:{1st_caller}_GT:{1st_caller}_DP:{1st_caller}_VD:{2nd_caller}_GT:{2nd_caller}_DP:{2nd_caller}_VD:{3rd_caller}_GT:{3rd_caller}_DP:{3rd_caller}_VD"
