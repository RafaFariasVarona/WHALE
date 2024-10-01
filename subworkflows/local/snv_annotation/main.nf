include { FORMAT2INFO            } from '../../../modules/local/format2info/main'
include { TABIX_TABIX            } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_BGZIP            } from '../../../modules/nf-core/tabix/bgzip/main'
include { AUTOMAP                } from '../../../modules/local/automap/main'
include { ENSEMBLVEP_VEP         } from '../../../modules/nf-core/ensemblvep/vep/main'

workflow SNV_ANNOTATION {

    take:
    merged_vcf
    fasta

    main:
    FORMAT2INFO (
        merged_vcf
    )

    TABIX_BGZIP (
        FORMAT2INFO.out.vcf_to_annotate
    )
    
    TABIX_TABIX (
        TABIX_BGZIP.out.output
    )

    sample_info = TABIX_BGZIP.out.output.join(TABIX_TABIX.out.tbi).join(FORMAT2INFO.out.fields)
    sample_info.view()

    AUTOMAP (
        merged_vcf,
		params.assembly,
		projectDir
    )

    AUTOMAP.out.roh_automap.view()

    ENSEMBLVEP_VEP (

        FORMAT2INFO.out.vcf_to_annotate,
        params.vep_genome, 
        params.vep_species, 
        params.vep_cache_version, 
        params.vep_cache, 
        fasta
        //params.vep_extra_files
    )


}
