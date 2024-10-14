#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import groovy.json.JsonSlurper
include { SANGERSEQ_BATCH } from './workflows/sangerseq_batch/main'
/*
 * pipeline input parameters
 */

log.info """\
    SANGER SEQUENCING DATA PROCESSING - P I P E L I N E
    ===================================================
    batch_id	: ${params.batch_id}
    trace_path  : ${params.trace_path}
    """
    .stripIndent()

workflow {
    //_2024-07-10-15-00-51_JO.ab1
    //2024_07_10_2024-07-10-15-00-51
    sample_trace_ch = Channel.fromFilePairs("${params.trace_path}/${params.batch_id}/*_{F,R}*_${params.trace_regex_suffix}", 
                checkIfExists:true)
    SANGERSEQ_BATCH ( sample_trace_ch )
}

workflow.onComplete {
    log.info ( workflow.success ? "\nSangerseq Analysis Workflow completed!\n" : "Oops .. something went wrong" )
}