#!/usr/bin/env nextflow

include { SANGERSEQ_ANALYSER  } from '../../modules/sangerseq/main'

workflow SANGERSEQ_BATCH {

    take:
    sample_trace_ch     // channel: [mandatory] [ val(sample), path(trace_files) ]

    main:
    ch_versions         = Channel.empty()

    SANGERSEQ_ANALYSER ( sample_trace_ch )
    ch_versions = ch_versions.mix(SANGERSEQ_ANALYSER.out.versions)

    emit:
    versions     = ch_versions                           // channel: [ path(versions.yml) ]
}