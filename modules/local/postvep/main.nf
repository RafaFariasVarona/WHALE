process POSTVEP {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(vep_tsv)
    tuple val(meta), path(roh_automap)
    path dbNSFP_gene
    path omim
    path regiondict
    path domino
    path tissue_expression
    val maf
    path genefilter
    path glowgenes
    val assembly
    path projectDir

    output:
    tuple val(meta), path("*.SNV.INDEL.annotated.tsv"), emit: pvm_tsv
    tuple val(meta), path("*.SNV.INDEL.annotated.tsv.xlsx"), emit: pvm_xlsx

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
	def omim_field = omim ? "--omim ${omim}" : ''
    def genefilter_field = genefilter ? "--genefilter ${genefilter}" : ''
	def glowgenes_field = glowgenes ? "--glowgenes ${glowgenes}" : ''

    """
    for vep_file in *vep.tsv
    do
        header_row="\$(head -n 1000 \${vep_file} | grep "#Uploaded_variation" -n | sed 's/:.*//')"

        Rscript ${projectDir}/bin/post-VEP_modification.R \\
        --input \${vep_file} \\
        --output \$(basename \${vep_file} .tsv).${assembly}.SNV.INDEL.annotated.tsv \\
        --numheader \${header_row} \\
        --dbNSFPgene ${dbNSFP_gene} \\
        --regiondict ${regiondict} \\
        --domino ${domino} \\
        --expression ${tissue_expression} \\
        --automap ./ \\
        --maf ${maf} \\
        ${omim_field} \\
        ${genefilter_field} \\
        ${glowgenes_field}
    done
    """
}
