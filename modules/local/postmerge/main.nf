process POSTMERGE {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(del_merged), path(inv_merged), path(ins_merged), path(bnd_merged), path(dup_merged)

    output:
    tuple val(meta), path("*merged_gt.bed"), emit: merged_gt
    tuple val(meta), path("*final.bed"), emit: merged_final

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    grep '#chr' \$(ls *.bed | head -n 1) > header.txt

    for sv_type in BND DEL DUP INS INV
    do
        cat *\${sv_type}*.bed | grep -v '#chr' | sort -k 1,1 -k2,2n | uniq > ${prefix}_\${sv_type}_tmp.bed
        cat header.txt ${prefix}_\${sv_type}_tmp.bed > ${prefix}_\${sv_type}_merged_gt_tmp.bed
        awk 'BEGIN{FS=OFS="\\t"}
        {
            for(i=1;i<=NF;i++)
                if(\$i=="")
                    \$i="."
        }
        NR==1 {
            \$26="sample_id"
        }
        NR>1 {
            \$26="${prefix}"
        } 1' ${prefix}_\${sv_type}_merged_gt_tmp.bed > ${prefix}_\${sv_type}_merged_gt.bed
    done

    grep '#chr' \$(ls *merged_gt.bed | head -n 1) > header.txt
    
    cat *merged_gt.bed | grep -v '#chr' | sort -k 1,1 -k2,2n > ${prefix}_merged_tmp.bed
    
    cat header.txt ${prefix}_merged_tmp.bed > ${prefix}_merged_final.bed
    """
}
