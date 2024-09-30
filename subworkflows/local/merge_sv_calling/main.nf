include { VCFTOBED1 as BEDSNIFFLES } from '../../../modules/local/vcftobed1/main'
include { VCFTOBED1 as BEDCUTESV   } from '../../../modules/local/vcftobed1/main'
include { VCFTOBED2 as BEDSVIM     } from '../../../modules/local/vcftobed2/main'
include { PREMERGE                 } from '../../../modules/local/premerge/main'
include { DEL_MERGE                } from '../../../modules/local/del_merge/main'
include { INV_MERGE                } from '../../../modules/local/inv_merge/main'
include { INS_MERGE                } from '../../../modules/local/ins_merge/main'
include { BND_MERGE                } from '../../../modules/local/bnd_merge/main'
include { DUP_MERGE                } from '../../../modules/local/dup_merge/main'
include { POSTMERGE                } from '../../../modules/local/postmerge/main'
include { ADD_SAMPLES_DB           } from '../../../modules/local/add_samples_db/main'

workflow MERGE_SV_CALLING {

    take:
    sniffles_vcf
    cutesv_vcf
    svim_vcf

    main:
    BEDSNIFFLES (
        sniffles_vcf
    )
    
    BEDCUTESV (
        cutesv_vcf
    )

    BEDSVIM (
        svim_vcf
    )
    
    bed_files = BEDSNIFFLES.out.bed.join(BEDCUTESV.out.bed).join(BEDSVIM.out.bed)
    
    PREMERGE (
        bed_files
    )
    
    DEL_MERGE (
        PREMERGE.out.sorted_del,
        PREMERGE.out.merged_del
    )
    
    INV_MERGE (
        PREMERGE.out.sorted_inv,
        PREMERGE.out.merged_inv
    )
    
    DUP_MERGE (
        PREMERGE.out.sorted_dup,
        PREMERGE.out.merged_dup
    )

    INS_MERGE (
        PREMERGE.out.sorted_ins
    )

    BND_MERGE (
        PREMERGE.out.sorted_bnd
    )

    merged_files = DEL_MERGE.out.del_merged_results.join(INV_MERGE.out.inv_merged_results).join(INS_MERGE.out.ins_merged_results).join(BND_MERGE.out.bnd_merged_results).join(DUP_MERGE.out.dup_merged_results)
    
    POSTMERGE (
        merged_files
    )

    if (params.sv_database == true) {
        
        multiinter = params.sv_multiinter ? Channel.fromPath(params.sv_multiinter).collect() : Channel.empty()

        ADD_SAMPLES_DB (
            POSTMERGE.out.merged_gt,
            multiinter
        )
    }

    emit:
    merged_final = params.sv_database == true ? ADD_SAMPLES_DB.out.samples_info : POSTMERGE.out.merged_final
}