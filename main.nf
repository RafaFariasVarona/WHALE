#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/variantpipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/variantpipeline
    Website: https://nf-co.re/variantpipeline
    Slack  : https://nfcore.slack.com/channels/variantpipeline
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { VARIANTPIPELINE  } from './workflows/variantpipeline'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_variantpipeline_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_variantpipeline_pipeline'

include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_variantpipeline_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO nf-core: Remove this line if you don't need a FASTA file
//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
// params.fasta = getGenomeAttribute('fasta')
// params.fasta_fai = getGenomeAttribute('fasta_fai')

// Initialize fasta file with meta map: from sarek/main.nf
fasta = params.fasta ? Channel.fromPath(params.fasta).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.empty()

// Initialize fasta_fai file with meta map:
fasta_fai = params.fasta_fai ? Channel.fromPath(params.fasta_fai).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.empty()

// Initialize fasta_gzi file with meta map:
fasta_gzi = params.fasta_gzi ? Channel.fromPath(params.fasta_gzi).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.empty()

// Initialize vep_extra_files

// vep plugins
vep_pluggin_files = []

vep_pluggin_files.add(file("${params.vep_annotation_dir}/${params.dbscSNV}", checkIfExists: true))
vep_pluggin_files.add(file("${params.vep_annotation_dir}/${params.dbscSNV_tbi}", checkIfExists: true))

vep_pluggin_files.add(file("${params.vep_annotation_dir}/${params.dbNSFP}", checkIfExists: true))
vep_pluggin_files.add(file("${params.vep_annotation_dir}/${params.dbNSFP_tbi}", checkIfExists: true))

vep_pluggin_files.add(file("${params.vep_annotation_gene_dir}/${params.loFtool}", checkIfExists: true))
vep_pluggin_files.add(file("${params.vep_annotation_gene_dir}/${params.exACpLI}", checkIfExists: true))

vep_pluggin_files.add(file("${params.vep_annotation_gene_dir}/${params.maxEntScan}", checkIfExists: true))

vep_pluggin_files.add(file("${params.vep_annotation_dir}/${params.cADD_INDELS}", checkIfExists: true))
vep_pluggin_files.add(file("${params.vep_annotation_dir}/${params.cADD_INDELS_tbi}", checkIfExists: true))

vep_pluggin_files.add(file("${params.vep_annotation_dir}/${params.cADD_SNVS}", checkIfExists: true))
vep_pluggin_files.add(file("${params.vep_annotation_dir}/${params.cADD_SNVS_tbi}", checkIfExists: true))

vep_pluggin_files_all = vep_pluggin_files ? Channel.fromPath(vep_pluggin_files).collect() : Channel.empty()
//vep_pluggin_files_all.view()

// vep custom

vep_custom_files = []

vep_custom_files.add(file("${params.vep_annotation_dir}/${params.kaviar}", checkIfExists: true))
vep_custom_files.add(file("${params.vep_annotation_dir}/${params.kaviar_tbi}", checkIfExists: true))

vep_custom_files.add(file("${params.vep_annotation_dir}/${params.mAF_FJD_COHORT}", checkIfExists: true))
vep_custom_files.add(file("${params.vep_annotation_dir}/${params.mAF_FJD_COHORT_tbi}", checkIfExists: true))

vep_custom_files.add(file("${params.vep_annotation_dir}/${params.cCRS_DB}", checkIfExists: true))
vep_custom_files.add(file("${params.vep_annotation_dir}/${params.cCRS_DB_tbi}", checkIfExists: true))

vep_custom_files.add(file("${params.vep_annotation_dir}/${params.cLINVAR}", checkIfExists: true))
vep_custom_files.add(file("${params.vep_annotation_dir}/${params.cLINVAR_tbi}", checkIfExists: true))

vep_custom_files.add(file("${params.vep_annotation_dir}/${params.dENOVO_DB}", checkIfExists: true))
vep_custom_files.add(file("${params.vep_annotation_dir}/${params.dENOVO_DB_tbi}", checkIfExists: true))

vep_custom_files.add(file("${params.vep_annotation_dir}/${params.gNOMADe_cov}", checkIfExists: true))
vep_custom_files.add(file("${params.vep_annotation_dir}/${params.gNOMADe_cov_tbi}", checkIfExists: true))
vep_custom_files.add(file("${params.vep_annotation_dir}/${params.gNOMADe}", checkIfExists: true))
vep_custom_files.add(file("${params.vep_annotation_dir}/${params.gNOMADe_tbi}", checkIfExists: true))
vep_custom_files.add(file("${params.vep_annotation_dir}/${params.gNOMADg_cov}", checkIfExists: true))
vep_custom_files.add(file("${params.vep_annotation_dir}/${params.gNOMADg_cov_tbi}", checkIfExists: true))
vep_custom_files.add(file("${params.vep_annotation_dir}/${params.gNOMADg}", checkIfExists: true))
vep_custom_files.add(file("${params.vep_annotation_dir}/${params.gNOMADg_tbi}", checkIfExists: true))

vep_custom_files.add(file("${params.vep_annotation_dir}/${params.mutScore}", checkIfExists: true))
vep_custom_files.add(file("${params.vep_annotation_dir}/${params.mutScore_tbi}", checkIfExists: true))

vep_custom_files.add(file("${params.vep_annotation_dir}/${params.spliceAI_SNV}", checkIfExists: true))
vep_custom_files.add(file("${params.vep_annotation_dir}/${params.spliceAI_SNV_tbi}", checkIfExists: true))

vep_custom_files.add(file("${params.vep_annotation_dir}/${params.spliceAI_INDEL}", checkIfExists: true))
vep_custom_files.add(file("${params.vep_annotation_dir}/${params.spliceAI_INDEL_tbi}", checkIfExists: true))

vep_custom_files.add(file("${params.vep_annotation_dir}/${params.REVEL}", checkIfExists: true))
vep_custom_files.add(file("${params.vep_annotation_dir}/${params.REVEL_tbi}", checkIfExists: true))

vep_custom_files.add(file("${params.vep_annotation_dir}/${params.CSVS_dir}/${params.cSVS}", checkIfExists: true))
vep_custom_files.add(file("${params.vep_annotation_dir}/${params.CSVS_dir}/${params.cSVS_tbi}", checkIfExists: true))

vep_custom_files_all = vep_custom_files ? Channel.fromPath(vep_custom_files).collect() : Channel.empty()
//vep_custom_files_all.view()



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_VARIANTPIPELINE {

    take:
    samplesheet // channel: samplesheet read in from --input
    fasta
    fasta_fai
    fasta_gzi
    vep_pluggin_files_all
    vep_custom_files_all
    
    main:

    //
    // WORKFLOW: Run pipeline
    //
    VARIANTPIPELINE (
        samplesheet,
        fasta,
        fasta_fai,
        fasta_gzi,
        vep_pluggin_files_all,
        vep_custom_files_all
    )

    emit:
    multiqc_report = VARIANTPIPELINE.out.multiqc_report // channel: /path/to/multiqc_report.html

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_VARIANTPIPELINE (
        PIPELINE_INITIALISATION.out.samplesheet,
        fasta,
        fasta_fai,
        fasta_gzi,
        vep_pluggin_files_all,
        vep_custom_files_all
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_VARIANTPIPELINE.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
