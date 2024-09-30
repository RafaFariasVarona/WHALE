include { CLAIR3                 } from '../../../modules/local/clair3/main'
include { FILTER_PASS            } from '../../../modules/local/filter/pass/main'
include { TABIX_TABIX            } from '../../../modules/nf-core/tabix/tabix/main'

workflow CLAIR3_SUB {

    take:
    ch_bambai
    fasta
    fasta_fai
    fasta_gzi

    main:
    CLAIR3 (
        ch_bambai,
        fasta,
        fasta_fai,
        fasta_gzi)
    
    FILTER_PASS (
        CLAIR3.out.vcf)

    TABIX_TABIX (
        FILTER_PASS.out.vcf)
    
    emit:
    vcf_tbi = FILTER_PASS.out.vcf.join(TABIX_TABIX.out.tbi) // channel: [ val(meta), path(vcf), path(tbi) ]
}