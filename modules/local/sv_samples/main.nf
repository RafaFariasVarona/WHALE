process SV_SAMPLES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    tuple val(meta), path(merged_sv)
    path(multiinter)

    output:
    tuple val(meta), path("*info_final.bed"), emit: samples_info
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    for sv_type in BND DEL DUP INS INV
    do
        bedtools map -a ${prefix}_\${sv_type}_merged_gt.bed -b \${sv_type}_multiinter.bed -c 5 -o collapse -header > ${prefix}_\${sv_type}_samples_info.bed
    done

    for file in *samples_info.bed
    do
        awk 'BEGIN{OFS="\\t"}
        NR==1 {print \$0, "overlapping_samples", "n_samples"; next} {
            samples = \$NF;
        
            for (i = 1; i < NF; i++) {
                printf "%s\\t", \$i;
            }

            split(samples, arr, ",");
            
            n = length(arr);
            for (i = 1; i <= n; i++) {
                for (j = i+1; j <= n; j++) {
                    if (arr[i] > arr[j]) {
                        temp = arr[i];
                        arr[i] = arr[j];
                        arr[j] = temp;
                    }
                }
            }

            unique = arr[1];
            for (i = 2; i <= n; i++) {
                if (arr[i] != arr[i-1]) {
                    unique = unique "," arr[i];
                }
            }

            n_samples = split(unique, tmp, ",");

            printf "%s\\t%d\\n", unique, n_samples;
        }' \${file} > \$(basename \${file} .bed)_tmp1.bed
    done

    grep '#chr' \$(ls *info_tmp1.bed | head -n 1) > header.txt

    cat ${prefix}*info_tmp1.bed | grep -v '#chr' | sort -k 1,1 -k2,2n > ${prefix}_samples_info_tmp2.bed
    
    cat header.txt ${prefix}_samples_info_tmp2.bed > ${prefix}_samples_info_final.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
