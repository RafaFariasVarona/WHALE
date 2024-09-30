include { DEEPVARIANT_SUB        } from '../deepvariant_sub/main'
include { NANOCALLER_SUB         } from '../nanocaller_sub/main'
include { CLAIR3_SUB             } from '../clair3_sub/main'

workflow SNV_CALLING {

    take:
    ch_bambai
    fasta
    fasta_fai
    fasta_gzi

    main:
    if (params.snv_caller == 'all') {
        DEEPVARIANT_SUB (
            ch_bambai,
            fasta,
            fasta_fai,
            fasta_gzi)
        
        NANOCALLER_SUB (
            ch_bambai,
            fasta,
            fasta_fai,
            fasta_gzi)
        
        CLAIR3_SUB (
            ch_bambai,
            fasta,
            fasta_fai,
            fasta_gzi)

    } else if (params.snv_caller == 'deepvariant') {
        DEEPVARIANT_SUB (
            ch_bambai,
            fasta,
            fasta_fai,
            fasta_gzi)

    } else if (params.snv_caller == 'nanocaller') {
        NANOCALLER_SUB (
            ch_bambai,
            fasta,
            fasta_fai,
            fasta_gzi)

    } else if (params.snv_caller == 'clair3') {
        CLAIR3_SUB (
            ch_bambai,
            fasta,
            fasta_fai,
            fasta_gzi)
    }
    
    emit:
    deepvariant_vcf_tbi = (params.snv_caller == 'all' || params.snv_caller == 'deepvariant') ? DEEPVARIANT_SUB.out.vcf_tbi : "" // channel: [ val(meta), path(vcf), path(tbi) ]
    nanocaller_vcf_tbi = (params.snv_caller == 'all' || params.snv_caller == 'nanocaller') ? NANOCALLER_SUB.out.vcf_tbi : "" // channel: [ val(meta), path(vcf), path(tbi) ]
    clair3_vcf_tbi = (params.snv_caller == 'all' || params.snv_caller == 'clair3') ? CLAIR3_SUB.out.vcf_tbi : "" // channel: [ val(meta), path(vcf), path(tbi) ]
}