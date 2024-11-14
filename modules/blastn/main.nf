// blastn
process BLASTN_QUERY {
    tag "$sample"
    label 'process_medium'

	input:
	tuple val(sample), path(contig_fasta)

	output:
    tuple val(sample), path("*_blast.asn")  , emit: asn
	tuple val(sample), path("*_blast.tsv")  , emit: tsv
	tuple val(sample), path("*_blast.html") , emit: html
    path "versions.yml"     				, emit: versions

	script:
	"""
	blastn \
		-query ${contig_fasta} -db ${params.blastn_db} \
        -remote \
		-outfmt 11 -out ${sample}_blast.asn \
		-evalue ${params.evalue} \
		-perc_identity ${params.perc_identity} \
		-qcov_hsp_perc ${params.qcov_hsp_perc} \
		-max_target_seqs ${params.max_target_seqs}

	blast_formatter \
		-archive ${sample}_blast.asn \
		-html -out ${sample}_blast.html

	blast_formatter \
		-archive ${sample}_blast.asn \
		-outfmt "7 saccver pident bitscore length evalue stitle" \
        -out ${sample}_blast.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        BLAST: \$(echo \$(blastn -version 2>&1))
    END_VERSIONS
	"""
}