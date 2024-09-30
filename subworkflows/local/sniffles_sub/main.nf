include { SNIFFLES               } from '../../../modules/nf-core/sniffles/main'
include { FILTER_PASS            } from '../../../modules/local/filter/pass/main'
include { BCFTOOLS_SORT          } from '../../../modules/nf-core/bcftools/sort/main'

workflow SNIFFLES_SUB {

    take:
    ch_bambai
    fasta
    fasta_fai
    fasta_gzi
    
    main:
    SNIFFLES (
        ch_bambai,
        fasta,
        fasta_fai,
        fasta_gzi)
    
    FILTER_PASS (
        SNIFFLES.out.vcf)
    
    BCFTOOLS_SORT (
        FILTER_PASS.out.vcf)
    
    emit:
    vcf = BCFTOOLS_SORT.out.vcf // channel: [ val(meta), path(vcf) ]
}
