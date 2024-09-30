include { CUTESV                 } from '../../../modules/nf-core/cutesv/main'
include { FILTER_PASS            } from '../../../modules/local/filter/pass/main'
include { BCFTOOLS_SORT          } from '../../../modules/nf-core/bcftools/sort/main'

workflow CUTESV_SUB {

    take:
    ch_bambai
    fasta
    fasta_fai
    fasta_gzi
    
    main:
    CUTESV (
        ch_bambai,
        fasta,
        fasta_fai,
        fasta_gzi)
    
    FILTER_PASS (
        CUTESV.out.vcf)
    
    BCFTOOLS_SORT (
        FILTER_PASS.out.vcf)
    
    emit:
    vcf = BCFTOOLS_SORT.out.vcf // channel: [ val(meta), path(vcf) ]
}
