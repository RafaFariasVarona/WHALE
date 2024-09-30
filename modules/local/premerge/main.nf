process PREMERGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    tuple val(meta), path(sniffles_bed), path(cutesv_bed), path(svim_bed)

    output:
    tuple val(meta), path("*het_DEL_merged.bed"), path("*hom_DEL_merged.bed"), emit: merged_del
    tuple val(meta), path("*het_INV_merged.bed"), path("*hom_INV_merged.bed"), emit: merged_inv
    tuple val(meta), path("*het_DUP_merged.bed"), path("*hom_DUP_merged.bed"), emit: merged_dup
    tuple val(meta), path("*het_DEL_sorted.bed"), path("*hom_DEL_sorted.bed"), emit: sorted_del
    tuple val(meta), path("*het_INV_sorted.bed"), path("*hom_INV_sorted.bed"), emit: sorted_inv
    tuple val(meta), path("*het_DUP_sorted.bed"), path("*hom_DUP_sorted.bed"), emit: sorted_dup
    tuple val(meta), path("*het_INS_sorted.bed"), path("*hom_INS_sorted.bed"), emit: sorted_ins
    tuple val(meta), path("*het_BND_sorted.bed"), path("*hom_BND_sorted.bed"), emit: sorted_bnd
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    for bed in ${sniffles_bed} ${cutesv_bed} ${svim_bed}
    do
        if grep -E -q '1/1|\\./\\.' \${bed}; then
            grep -E '1/1|\\./\\.' \${bed} > \$(basename \${bed} _filtered_sorted.bed)_hom.bed
        else
            echo "1/1 genotype not found in \${bed}"
        fi
    
        if grep -E -q '0/1|0/0|\\./\\.' \${bed}; then
            grep -E '0/1|0/0|\\./\\.' \${bed} > \$(basename \${bed} _filtered_sorted.bed)_het.bed
        else
            echo "0/1 or 0/0 genotype not found in \${bed}"
        fi
    done

    for file in *het.bed *hom.bed
    do
        for sv_type in BND DEL DUP INS INV DUP_TANDEM DUP_INT
        do
            if grep -q "\\.\${sv_type}\\." \${file}; then
                grep "\\.\${sv_type}\\." \${file} > \$(basename \${file} .bed)_\${sv_type}.bed
            else
                echo "SV type \${sv_type} not found in \${file}"
            fi
        done
    done

    for GT in hom het
    do
        for sv_type in DEL INV INS BND DUP
        do
            file_count=\$(ls ${prefix}*\${GT}_\${sv_type}*.bed 2>/dev/null | wc -l)
            if [ "\${file_count}" -ge 2 ]; then
                cat ${prefix}*\${GT}_\${sv_type}*.bed > ${prefix}_\${GT}_\${sv_type}_pasted.bed
            elif [ "\${file_count}" -eq 1 ]; then
                mv ${prefix}*\${GT}_\${sv_type}*.bed ${prefix}_\${GT}_\${sv_type}_pasted.bed
            fi
            if [ -f "${prefix}_\${GT}_\${sv_type}_pasted.bed" ]; then
                sort -k 1,1 -k2,2n ${prefix}_\${GT}_\${sv_type}_pasted.bed > ${prefix}_\${GT}_\${sv_type}_sorted.bed
                if [[ "\${sv_type}" == "DEL" || "\${sv_type}" == "INV" || "\${sv_type}" == "DUP" ]]; then
                    bedtools merge -c 4 -o collapse -i ${prefix}_\${GT}_\${sv_type}_sorted.bed > ${prefix}_\${GT}_\${sv_type}_merged.bed
                fi
            fi
        done
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
