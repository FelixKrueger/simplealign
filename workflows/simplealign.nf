/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowSimplealign.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { MULTIQC                     }      from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS }      from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE } from '../subworkflows/nf-core/fastq_fastqc_umitools_trimgalore/main'
include { FASTQ_ALIGN_BOWTIE2 }              from '../subworkflows/nf-core/fastq_align_bowtie2/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow SIMPLEALIGN {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input)
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema


    //
    // SUBWORKFLOW: Read QC and trim adapters
    //

    // take:
    // reads             // channel: [ val(meta), [ reads ] ]
    // skip_fastqc       // boolean: true/false
    // with_umi          // boolean: true/false
    // skip_umi_extract  // boolean: true/false
    // skip_trimming     // boolean: true/false
    // umi_discard_read  // integer: 0, 1 or 2
    // min_trimmed_reads // integer: > 0

    FASTQ_FASTQC_UMITOOLS_TRIMGALORE (
        INPUT_CHECK.out.reads,
        params.skip_fastqc,
        false,                      // with_umi          // boolean: true/false
        true,                       // skip_umi_extract  // boolean: true/false
        params.skip_trimming,
        1,                          // umi_discard_read  // integer: 0, 1 or 2
        1000                        // min_trimmed_reads // integer: > 0
    ) 
    ch_filtered_reads      = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.reads
    ch_fastqc_raw_multiqc  = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip
    ch_fastqc_trim_multiqc = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip
    ch_trim_log_multiqc    = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_log
    ch_trim_read_count     = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_read_count
    ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.versions)


    //
    // SUBWORKFLOW: Alignment with Bowtie2 & BAM QC
    //

    // This is what the subworkflow takes as arguments:

    // workflow FASTQ_ALIGN_BOWTIE2 {
    // take:
    // ch_reads          // channel: [ val(meta), [ reads ] ]
    // ch_index          // channel: /path/to/bowtie2/index/
    // save_unaligned    // val
    // sort_bam          // val
    // ch_fasta          // channel: /path/to/reference.fasta
    
    //FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.reads.view()
    
    FASTQ_ALIGN_BOWTIE2 (
        FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.reads,
        params.genomes[params.genome]['bowtie2'], // assuming the index has been built already
        params.save_unaligned,
        false,
        params.genomes[params.genome]['fasta']
    )
    ch_genome_bam        = FASTQ_ALIGN_BOWTIE2.out.bam
    ch_genome_bam_index  = FASTQ_ALIGN_BOWTIE2.out.bai
    ch_samtools_stats    = FASTQ_ALIGN_BOWTIE2.out.stats
    ch_samtools_flagstat = FASTQ_ALIGN_BOWTIE2.out.flagstat
    ch_samtools_idxstats = FASTQ_ALIGN_BOWTIE2.out.idxstats
    ch_versions = ch_versions.mix(FASTQ_ALIGN_BOWTIE2.out.versions.first())
    
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )


    //
    // MODULE: MultiQC
    //

    // input:
    // path  multiqc_files, stageAs: "?/*"
    // path(multiqc_config)
    // path(extra_multiqc_config)
    // path(multiqc_logo)

    workflow_summary    = WorkflowSimplealign.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)
    
    methods_description    = WorkflowSimplealign.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)
    
    // ch_multiqc_files = Channel.empty()
    // ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]))
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]))
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]))
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BOWTIE2.out.log.collect{it[1]}.ifEmpty([]))
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BOWTIE2.out.zip.collect{it[1]}.ifEmpty([]))
    // ch_multiqc_files.view()   
    ch_fastqc_raw_multiqc.collect{it[1]}.ifEmpty([]).view()

    // ::::::::::::::::: FROM SEQC :::::::::::::::::::
    workflow_summary    = WorkflowSeqc.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowSeqc.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    
    // ::::::::::::::::: FROM SEQC ::::::::::::::::::

    // MULTIQC (
    //     ch_multiqc_config,
    //     ch_multiqc_custom_config.collect().ifEmpty([]),
    //     CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect(),
    //     ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
    //     ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
    //     ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'),
    //     ch_multiqc_logo.collect().ifEmpty([]),
    //     ch_fastqc_raw_multiqc.collect{it[1]}.ifEmpty([]),
    //     ch_fastqc_trim_multiqc.collect{it[1]}.ifEmpty([]),
    //     ch_trim_log_multiqc.collect{it[1]}.ifEmpty([]),
    //     ch_samtools_stats.collect{it[1]}.ifEmpty([]),
    //     ch_samtools_flagstat.collect{it[1]}.ifEmpty([]),
    //     ch_samtools_idxstats.collect{it[1]}.ifEmpty([])
    //     // ch_multiqc_files.collect(),
    //     // ch_multiqc_config.toList(),
    //     // ch_multiqc_custom_config.toList().ifEmpty([]),
    //     // ch_multiqc_logo.collect().ifEmpty([]),
    // )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
