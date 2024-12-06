#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//VALIDATE INPUTS
def checkParamList = [
    params.output_dir,
    params.trace_path,
    params.batch_id,
    params.trace_regex_suffix]

for (param in checkParamList) if (!param) error("Required options were not provided")

include { SANGERSEQ_BATCH } from './workflows/sangerseq_batch/main'
include { BLASTN_SEARCH   } from './workflows/blastn_search/main'
/*
 * pipeline input parameters
 */

log.info """\
    SANGER SEQUENCING DATA PROCESSING - P I P E L I N E
    ===================================================
    batch_id            : ${params.batch_id}
    trace_path          : ${params.trace_path}
    trace_regex_suffix  : ${params.trace_regex_suffix}
    trim_cutoff         : ${params.trim_cutoff}
    min_seq_len         : ${params.min_seq_len}
    output_dir          : ${params.output_dir}
    perc_identity       : ${params.perc_identity}
    max_target_seqs     : ${params.max_target_seqs}
    """
    .stripIndent()

workflow {
    ch_versions = Channel.empty()
    sample_trace_ch = Channel.fromFilePairs("${params.trace_path}/${params.batch_id}/*_{F,R}*_${params.trace_regex_suffix}", 
                checkIfExists:true)
    SANGERSEQ_BATCH ( sample_trace_ch )
    ch_versions = ch_versions.mix( SANGERSEQ_BATCH.out.versions )
    BLASTN_SEARCH(SANGERSEQ_BATCH.out.contig_fasta)
    ch_versions = ch_versions.mix( BLASTN_SEARCH.out.versions )
}

workflow.onComplete {
    log.info ( workflow.success ? "\nSanger sequencing batch processing completed.\n" : "Oops .. something went wrong" )
}