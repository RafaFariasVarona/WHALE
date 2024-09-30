process NANOCALLER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanocaller:3.5.0--hdfd78af_0':
        'biocontainers/nanocaller:3.5.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(gzi)

    output:
    tuple val(meta), path("*nanocaller.vcf.gz"), emit: vcf
    tuple val(meta), path("*nanocaller.vcf.gz.tbi"), emit: tbi

    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    NanoCaller \\
        --bam ${bam} \\
        --ref ${fasta} \\
        --cpu ${task.cpus} \\
        --prefix ${prefix}_nanocaller \\
        --sample ${prefix}.NC \\
        ${args}
    
    mv ${prefix}_nanocaller.vcf.gz ${prefix}_nanocaller_tmp.vcf.gz
    
    gzip -dc ${prefix}_nanocaller_tmp.vcf.gz | sed "s/Number=\\./Number=R/g" > ${prefix}_nanocaller_tmp.vcf
        
    grep -E -v 'ORP|ORL|GQ' ${prefix}_nanocaller_tmp.vcf > ${prefix}_nanocaller.vcf

    gzip ${prefix}_nanocaller.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        NanoCaller: \$(echo \$(NanoCaller --help) | grep Version)
    END_VERSIONS
    """
}
// Para hacer el "merge" de los vcfs, es necesario que en el header de los vcfs ponga ##FORMAT=<ID=AD,Number=R
// En el vcf de NanoCaller, en AD pone "Number=.", por lo que bcftools merge no funciona. Para que funcione, se
// inclue el comando sed para cambiar "Number=." por "Number=R".
// En el vcf de NanoCaller, aparecen unas "variantes" definidas por ORP, ORL y GQ que dan problemas en el merge.
// Para filtrar esas "variantes" se incluye el comando grep -v