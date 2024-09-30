process SOFTWARE_INFO {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(GT)
    tuple val(meta), path(caller_order)

    output:
    tuple val(meta), path("${prefix}_SF.txt"), emit: sf

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    cut -f 1 ${GT} | sed 's/.+/'"\$(awk 'NR==1 {print}' ${caller_order})"'/' -r > ${prefix}_1.txt
    cut -f 2 ${GT} | sed 's/.+/'"\$(awk 'NR==2 {print}' ${caller_order})"'/' -r > ${prefix}_2.txt
    cut -f 3 ${GT} | sed 's/.+/'"\$(awk 'NR==3 {print}' ${caller_order})"'/' -r > ${prefix}_3.txt

    paste -d "\\t" ${prefix}_1.txt ${prefix}_2.txt ${prefix}_3.txt | perl -alne 'print join "_", @F' > ${prefix}_SF.txt
    """
}
// Get software information (DV, C3, NC) according to presence or absence of genotype