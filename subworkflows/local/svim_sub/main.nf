include { SVIM                   } from '../../../modules/local/svim/main'
include { FILTER_QUAL            } from '../../../modules/local/filter/qual/main'
include { BCFTOOLS_SORT          } from '../../../modules/nf-core/bcftools/sort/main'

workflow SVIM_SUB {

    take:
    ch_bambai
    fasta
    fasta_fai
    fasta_gzi
    
    main:
    SVIM (
        ch_bambai,
        fasta,
        fasta_fai,
        fasta_gzi)
    
    FILTER_QUAL (
        SVIM.out.vcf)
    
    BCFTOOLS_SORT (
        FILTER_QUAL.out.vcf)
    
    emit:
    vcf = BCFTOOLS_SORT.out.vcf // channel: [ val(meta), path(vcf) ]
}
