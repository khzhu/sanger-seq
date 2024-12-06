#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { BLASTN_QUERY  } from '../../modules/blastn/main'
include { BLASTN_HTML   } from '../../modules/blast_html/main'
include { SIGNOUT_FORM  } from '../../modules/signout_form/main'

workflow BLASTN_SEARCH {

    take:
    contig_fasta_ch // channel: [mandatory] [ val(sample), path(fasta_file) ]

    main:
    ch_versions = Channel.empty()

    BLASTN_QUERY ( contig_fasta_ch )
    ch_versions = ch_versions.mix(BLASTN_QUERY.out.versions)
    BLASTN_HTML ( BLASTN_QUERY.out.tsv, BLASTN_QUERY.out.html )
    ch_versions = ch_versions.mix(BLASTN_HTML.out.versions)
    SIGNOUT_FORM ( BLASTN_HTML.out.html )

    emit:
    asn      = BLASTN_QUERY.out.asn
    tsv      = BLASTN_QUERY.out.tsv
    html     = BLASTN_HTML.out.html
    docx     = SIGNOUT_FORM.out.docx
    versions = ch_versions  // channel: [ path(versions.yml) ]
}