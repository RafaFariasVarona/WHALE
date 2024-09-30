include { ANNOTSV_INSTALLANNOTATIONS } from '../../../modules/nf-core/annotsv/installannotations/main'
include { ANNOTSV_ANNOTSV            } from '../../../modules/nf-core/annotsv/annotsv/main' 

workflow SV_ANNOTATION {

    take:
    merged_bed

    main:
    if (params.annotsv_annotations == 'install') {
        ANNOTSV_INSTALLANNOTATIONS()
        
        annotations = ANNOTSV_INSTALLANNOTATIONS.out.annotations
    } else {
        annotations = params.annotsv_annotations ? Channel.fromPath(params.annotsv_annotations).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.empty()
    }

    gene_transcripts = params.gene_transcripts ? Channel.fromPath(params.gene_transcripts).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.empty()

    ANNOTSV_ANNOTSV (
        merged_bed,
        annotations,
        gene_transcripts
    )

    emit:
    annotated_tsv = ANNOTSV_ANNOTSV.out.tsv
}
