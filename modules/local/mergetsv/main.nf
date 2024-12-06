process MERGETSV {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(tsv)
    val assembly

    output:
    tuple val(meta), path("*final.tsv"), emit: final_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    head -n 1 *split00*.tsv > header.txt

    tail -q -n +2 *.tsv | cat > merged.tsv

    cat header.txt merged.tsv > ${prefix}.${assembly}.SNV.INDEL.annotated.final.tsv
    """
}
