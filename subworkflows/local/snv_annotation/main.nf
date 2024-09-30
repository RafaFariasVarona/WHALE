include { FORMAT2INFO            } from '../../../modules/local/format2info/main'
include { TABIX_TABIX            } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_BGZIP            } from '../../../modules/nf-core/tabix/bgzip/main'
include { AUTOMAP                } from '../../../modules/local/automap/main'

workflow SNV_ANNOTATION {

    take:
    merged_vcf

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
}
