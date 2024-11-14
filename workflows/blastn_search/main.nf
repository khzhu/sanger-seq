#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { BLASTN_QUERY  } from '../../modules/blastn/main'

workflow BLASTN_SEARCH {

    take:
    contig_fasta_ch             // channel: [mandatory] [ val(sample), path(fasta_file) ]

    main:
    ch_versions         = Channel.empty()

    BLASTN_QUERY ( contig_fasta_ch )
    ch_versions = ch_versions.mix(BLASTN_QUERY.out.versions)

    emit:
    asn         = BLASTN_QUERY.out.asn
    tsv         = BLASTN_QUERY.out.tsv
    html        = BLASTN_QUERY.out.html
    versions    = ch_versions  // channel: [ path(versions.yml) ]
}