include { NANOCALLER             } from '../../../modules/local/nanocaller/main'
include { FILTER_PASS            } from '../../../modules/local/filter/pass/main'
include { TABIX_TABIX            } from '../../../modules/nf-core/tabix/tabix/main'

workflow NANOCALLER_SUB {

    take:
    ch_bambai
    fasta
    fasta_fai
    fasta_gzi

    main:
    NANOCALLER (
        ch_bambai,
        fasta,
        fasta_fai,
        fasta_gzi)
    
    FILTER_PASS (
        NANOCALLER.out.vcf)
    
    TABIX_TABIX (
        FILTER_PASS.out.vcf)

    emit:
    vcf_tbi = FILTER_PASS.out.vcf.join(TABIX_TABIX.out.tbi) // channel: [ val(meta), path(vcf), path(tbi) ]
}