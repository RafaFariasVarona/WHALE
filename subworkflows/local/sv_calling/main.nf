include { SNIFFLES_SUB           } from '../sniffles_sub/main'
include { CUTESV_SUB             } from '../cutesv_sub/main'
include { SVIM_SUB               } from '../svim_sub/main'

workflow SV_CALLING {

    take:
    ch_bambai
    fasta
    fasta_fai
    fasta_gzi
    
    main:
    if (params.sv_caller == 'all') {
        SNIFFLES_SUB (
            ch_bambai,
            fasta,
            fasta_fai,
            fasta_gzi)

        CUTESV_SUB (
            ch_bambai,
            fasta,
            fasta_fai,
            fasta_gzi)
        
        SVIM_SUB (
            ch_bambai,
            fasta,
            fasta_fai,
            fasta_gzi) 
        
    }
    else if (params.sv_caller == 'sniffles') {
        SNIFFLES_SUB (
            ch_bambai,
            fasta,
            fasta_fai,
            fasta_gzi)

    }
    else if (params.sv_caller == 'cutesv') {
        CUTESV_SUB (
            ch_bambai,
            fasta,
            fasta_fai,
            fasta_gzi)
        
    }    
    else if (params.sv_caller == 'svim') {
        SVIM_SUB (
            ch_bambai,
            fasta,
            fasta_fai,
            fasta_gzi)

    }
    emit:
    sniffles_vcf = (params.sv_caller == 'all' || params.sv_caller == 'sniffles') ? SNIFFLES_SUB.out.vcf : "" // channel: [ val(meta), path(vcf) ]
    cutesv_vcf = (params.sv_caller == 'all' || params.sv_caller == 'cutesv') ? CUTESV_SUB.out.vcf : "" // channel: [ val(meta), path(vcf) ]
    svim_vcf = (params.sv_caller == 'all' || params.sv_caller == 'svim') ? SVIM_SUB.out.vcf : "" // channel: [ val(meta), path(vcf) ]
}
